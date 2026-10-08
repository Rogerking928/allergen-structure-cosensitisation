"""把相關性推到判別力：結構相似度與序列相似度各能認出多少「強共同致敏」的成分對。

為什麼要做這件事：相關係數（Spearman ρ）回答「排序一致嗎」，
但篩選準則（Codex 的 80 殘基窗 35%）實際上是個二分判定，審稿人要看的是
敏感度／特異度與可否改善，不是 ρ。本檔算三件事：

1. 判別力（AUC）：以 z > 3 當「強共同致敏」，比較 TM-score、80 殘基窗 identity、
   全長 identity。**這不是外部驗證**——z 來自同一批估計，AUC 只描述兩個相似度
   在這份資料裡認出強訊號的能力，寫稿時必須講明。
2. 兩個門檻的 2×2：Codex 序列準則 vs 結構門檻 TM ≥ 0.5，各抓到多少、重疊多少。
3. 落在「序列不達標、結構達標」那一格的對逐一列出（具名，供臨床對照）。

信賴區間一律用**節點自助法**（resample 過敏原成分，不是 resample 對）：
同一個成分出現在很多對裡，對不獨立，直接 bootstrap 對會把區間算得太窄。
"""
import json
import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score, roc_curve
from common import ROOT
import structure_cosens as S

OUT = ROOT / "out"
RNG = np.random.default_rng(20261007)
Z_HI = 3.0          # 「強共同致敏」：合併後 z > 3
N_BOOT = 2000


def core_pairs():
    d = pd.read_csv(OUT / "pairs_structure_cosens.csv")
    c = d[(~d.same_species) & (~d.subchain) & (d.len_a >= S.MIN_LEN) & (d.len_b >= S.MIN_LEN)
          & d.annotated].copy()
    c["strong"] = (c.z > Z_HI).astype(int)
    c["codex"] = c.seqid80 >= S.SEQ_HI
    c["struct"] = c.tm >= S.TM_HI
    return c


def node_boot_weights(g, rng):
    """節點自助法：成分抽後放回，一對的權重＝兩端被抽中的次數相乘。"""
    nodes = pd.unique(pd.concat([g.a, g.b]))
    cnt = pd.Series(rng.multinomial(len(nodes), np.full(len(nodes), 1 / len(nodes))), index=nodes)
    return cnt.reindex(g.a).values * cnt.reindex(g.b).values


def auc_ci(g, col, n=N_BOOT):
    y = g.strong.values
    x = g[col].values
    obs = roc_auc_score(y, x)
    rng = np.random.default_rng(RNG.integers(1 << 32))
    vals = []
    for _ in range(n):
        w = node_boot_weights(g, rng)
        if w.sum() == 0:
            continue
        yy = y[w > 0]
        if yy.min() == yy.max():
            continue
        vals.append(roc_auc_score(y, x, sample_weight=w))
    lo, hi = np.percentile(vals, [2.5, 97.5])
    return dict(predictor=col, n=len(g), events=int(y.sum()), auc=obs, lo=lo, hi=hi)


def point_ci(g, mask_col, n=N_BOOT):
    """二分準則的敏感度／特異度，同樣用節點自助法。"""
    y = g.strong.values.astype(bool)
    m = g[mask_col].values.astype(bool)
    sens, spec = (m & y).sum() / y.sum(), (~m & ~y).sum() / (~y).sum()
    rng = np.random.default_rng(RNG.integers(1 << 32))
    se, sp = [], []
    for _ in range(n):
        w = node_boot_weights(g, rng)
        if (w * y).sum() == 0 or (w * ~y).sum() == 0:
            continue
        se.append((w * m * y).sum() / (w * y).sum())
        sp.append((w * ~m * ~y).sum() / (w * ~y).sum())
    return dict(criterion=mask_col, flagged=int(m.sum()),
                sens=sens, sens_lo=np.percentile(se, 2.5), sens_hi=np.percentile(se, 97.5),
                spec=spec, spec_lo=np.percentile(sp, 2.5), spec_hi=np.percentile(sp, 97.5))


def main():
    c = core_pairs()
    hom = c[c.homologous]
    res = {"z_threshold": Z_HI, "n_core": len(c), "n_strong": int(c.strong.sum()),
           "n_homologous": len(hom), "n_strong_homologous": int(hom.strong.sum())}

    res["auc_all"] = [auc_ci(c, col) for col in ("tm", "seqid80", "seqid")]
    res["auc_homologous"] = [auc_ci(hom, col) for col in ("tm", "seqid80")]
    res["points_all"] = [point_ci(c, m) for m in ("codex", "struct")]

    # 2×2：兩個門檻怎麼分配這 8,239 對
    tab = []
    for st in (True, False):
        for cx in (True, False):
            g = c[(c.struct == st) & (c.codex == cx)]
            tab.append(dict(structural=st, codex=cx, n=len(g), strong=int(g.strong.sum()),
                            pct_strong=100 * g.strong.mean() if len(g) else np.nan,
                            median_or=float(np.exp(g.beta.median())) if len(g) else np.nan))
    res["reclassification"] = tab
    pd.DataFrame(tab).to_csv(OUT / "reclassification.csv", index=False)

    # 關鍵格：序列不達標但結構達標
    q = c[c.struct & ~c.codex].sort_values("z", ascending=False)
    q[["a", "b", "tm", "seqid80", "beta", "se", "z", "fam_a", "fam_b", "org_a", "org_b",
       "homologous"]].to_csv(OUT / "quadrant_struct_only.csv", index=False)
    res["quadrant_struct_only"] = dict(
        n=len(q), strong=int(q.strong.sum()), median_or=float(np.exp(q.beta.median())),
        pct_homologous=100 * q.homologous.mean())
    bg = c[~c.struct & ~c.codex]
    res["quadrant_neither"] = dict(n=len(bg), strong=int(bg.strong.sum()),
                                   median_or=float(np.exp(bg.beta.median())))

    # ROC 曲線座標（給圖用）
    for col in ("tm", "seqid80"):
        fpr, tpr, _ = roc_curve(c.strong, c[col])
        pd.DataFrame({"fpr": fpr, "tpr": tpr}).to_csv(OUT / f"roc_{col}.csv", index=False)

    (OUT / "discrimination.json").write_text(json.dumps(res, indent=1))
    for r in res["auc_all"]:
        print(f"all pairs  {r['predictor']:>8}  AUC {r['auc']:.3f} ({r['lo']:.3f}-{r['hi']:.3f})")
    for r in res["auc_homologous"]:
        print(f"homologous {r['predictor']:>8}  AUC {r['auc']:.3f} ({r['lo']:.3f}-{r['hi']:.3f})")
    for r in res["points_all"]:
        print(f"{r['criterion']:>8}: flags {r['flagged']}, sens {r['sens']:.2f} "
              f"({r['sens_lo']:.2f}-{r['sens_hi']:.2f}), spec {r['spec']:.3f}")
    print(pd.DataFrame(tab).to_string(index=False))


if __name__ == "__main__":
    main()
