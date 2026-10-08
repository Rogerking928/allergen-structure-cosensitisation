"""把「強共同致敏」改在**另一半病人**身上定義，拿掉結局與預測子共用同一批估計的循環。

主分析的判別力是同一個世代內的一致性：z 來自那批勝算比，TM 又拿去預測它。
TM-score 完全不用病人資料，所以只要把**結局**換成一組不相交的病人，循環就消失。

做法：在每個晶片平台內把病人隨機對半分（分層以免某一半缺整個平台），
兩半各自獨立重算全部成分對的校正後勝算比，然後
  1. 兩半各自的 ρ(TM, log OR)——穩定性；
  2. 以 B 半定義的強共同致敏（|z| > 3）當結局、TM 與序列當預測子算 AUC——**樣本外**；
  3. 兩半 log OR 的一致性（Spearman 與符號一致率）。
A、B 兩個方向都做，重複 R 次取中位數與範圍。
"""
import itertools
import json
import sys
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.metrics import roc_auc_score
from common import load, ROOT
import cosens
import structure_cosens as S

OUT = ROOT / "out"
R = 5
SEED = 20261007


def estimate(d, comps):
    """cosens.main 的核心，抽出來給半樣本用（同樣的閾值與合併方式）。"""
    out = []
    for chip in d.Chip_Type.unique():
        y, keep = cosens.strat_table(d, comps, chip)
        if len(y) < 100:
            continue
        tot = y.sum(axis=1).values
        npos = y.sum(axis=0)
        use = [c for c in keep if npos[c] >= cosens.MIN_POS]
        arr = y[use].values
        idx = {c: k for k, c in enumerate(use)}
        for a, b in itertools.combinations(use, 2):
            ya, yb = arr[:, idx[a]], arr[:, idx[b]]
            m = ~(np.isnan(ya) | np.isnan(yb))
            if m.sum() < 50 or ((ya[m] == 1) & (yb[m] == 1)).sum() < cosens.MIN_CELL:
                continue
            breadth = tot[m] - ya[m] - yb[m]
            beta, se, how = cosens.fit_pair(ya[m], yb[m], breadth)
            if how != "mle" or not np.isfinite(se):
                continue
            out.append(dict(a=a, b=b, beta=beta, se=se, w=1 / se ** 2))
    r = pd.DataFrame(out)
    g = r.groupby(["a", "b"])
    p = pd.DataFrame({"beta": g.apply(lambda x: np.average(x.beta, weights=x.w),
                                      include_groups=False),
                      "se": 1 / np.sqrt(g.w.sum())}).reset_index()
    p["z"] = p.beta / p.se
    return p


def core_keys():
    """主分析那 8,239 對（同源與否、TM、序列相似度都來自結構側，與病人無關）。"""
    d = pd.read_csv(OUT / "pairs_structure_cosens.csv")
    return d[(~d.same_species) & (~d.subchain) & (d.len_a >= S.MIN_LEN)
             & (d.len_b >= S.MIN_LEN) & d.annotated][
        ["a", "b", "tm", "seqid80", "homologous"]]


def split(d, rng):
    half = np.zeros(len(d), dtype=int)
    for chip in d.Chip_Type.unique():
        i = np.where(d.Chip_Type.values == chip)[0]
        perm = rng.permutation(i)
        half[perm[:len(i) // 2]] = 1
    return half


def main():
    reps = int(sys.argv[1]) if len(sys.argv) > 1 else R
    d, meta, comps = load()
    key = core_keys()
    rng = np.random.default_rng(SEED)
    rows = []
    for rep in range(reps):
        h = split(d, rng)
        est = [estimate(d[h == k], comps) for k in (0, 1)]
        both = est[0].merge(est[1], on=["a", "b"], suffixes=("_0", "_1"))
        m = key.merge(both, on=["a", "b"])
        hom = m[m.homologous]
        agree = dict(rep=rep, n_pairs=len(m), n_homologous=len(hom),
                     rho_halves=stats.spearmanr(m.beta_0, m.beta_1).statistic,
                     sign_agreement=float((np.sign(m.beta_0) == np.sign(m.beta_1)).mean()))
        for k in (0, 1):
            agree[f"rho_tm_half{k}"] = stats.spearmanr(hom.tm, hom[f"beta_{k}"]).statistic
        # 樣本外：結局在另一半定義
        for k, other in ((0, 1), (1, 0)):
            y = (m[f"z_{other}"] > 3).astype(int)   # 與主分析一致：單側，不取絕對值
            if y.min() == y.max():
                continue
            agree[f"auc_tm_out{other}"] = roc_auc_score(y, m.tm)
            agree[f"auc_seq_out{other}"] = roc_auc_score(y, m.seqid80)
            flag = m.tm >= S.TM_HI
            cx = m.seqid80 >= S.SEQ_HI
            agree[f"struct_only_n_out{other}"] = int((flag & ~cx).sum())
            agree[f"struct_only_strong_out{other}"] = int((flag & ~cx & (y == 1)).sum())
            agree[f"neither_pct_out{other}"] = float(100 * y[~flag & ~cx].mean())
        rows.append(agree)
        print(f"rep {rep}: pairs {len(m)}, rho halves {agree['rho_halves']:.2f}, "
              f"rho TM {agree['rho_tm_half0']:.2f}/{agree['rho_tm_half1']:.2f}, "
              f"AUC out {agree.get('auc_tm_out1', float('nan')):.2f}", flush=True)
    t = pd.DataFrame(rows)
    t.to_csv(OUT / "holdout.csv", index=False)

    def rng_(col):
        v = pd.concat([t[c] for c in t.columns if c.startswith(col)])
        return dict(median=float(v.median()), lo=float(v.min()), hi=float(v.max()))

    res = dict(reps=reps, n_pairs=int(t.n_pairs.median()),
               n_homologous=int(t.n_homologous.median()),
               rho_between_halves=rng_("rho_halves"),
               sign_agreement=rng_("sign_agreement"),
               rho_tm=rng_("rho_tm_half"),
               auc_tm_out_of_sample=rng_("auc_tm_out"),
               auc_seq_out_of_sample=rng_("auc_seq_out"),
               struct_only_n=rng_("struct_only_n_out"),
               struct_only_strong=rng_("struct_only_strong_out"),
               neither_pct=rng_("neither_pct_out"))
    (OUT / "holdout.json").write_text(json.dumps(res, indent=1))
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
