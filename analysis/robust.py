"""同源對內「結構相似度 → 共同致敏」的穩健性檢查。

1. 結構是否在序列之外還有增量：兩者先轉秩，log OR 對 TM 與序列 identity 的秩做迴歸，
   顯著性用過敏原標籤排列（同 structure_cosens.mantel 的理由：對不是獨立的）。
2. 留一家族：每次拿掉一個同源群，看 Spearman 是否仍成立（排除單一家族撐起全部結果）。
3. 分晶片：ISAC v1／ISAC v2／ALEX 各自的 log OR 分別重算。
4. 門檻：TM 0.4／0.5／0.6／0.7 各區間的 log OR 中位數（給臨床讀者的對照表）。
"""
import numpy as np
import pandas as pd
from scipy import stats
from common import ROOT
import homology
import structure_cosens as S

OUT = ROOT / "out"
RNG = np.random.default_rng(20261003)
# 序列相似度一律用 80-aa 窗 identity（Codex 準則、Figure 3、AUC 用的同一個量）。
# 2026-10-08 以前這裡用全長 identity（seqid），與稿件其他地方不一致，被審閱抓到
SEQ = "seqid80"
# 哺乳類的目名不照 -ales／-formes 結尾，逐一列出
MAMMAL_ORDERS = {"Carnivora", "Rodentia", "Primates", "Artiodactyla", "Perissodactyla",
                 "Lagomorpha", "Eulipotyphla"}


def homologous_core():
    d = pd.read_csv(OUT / "pairs_structure_cosens.csv")
    core = d[(~d.same_species) & (~d.subchain) & (d.len_a >= S.MIN_LEN) & (d.len_b >= S.MIN_LEN)
             & d.annotated & d.homologous].copy()
    return core


def incremental(g, n=5000):
    """log OR ~ rank(TM) + rank(80-aa identity)；以過敏原標籤排列 log OR 求 TM 係數的 p。"""
    def coef(y, x1, x2):
        X = np.column_stack([np.ones(len(y)), stats.rankdata(x1), stats.rankdata(x2)])
        X[:, 1:] = (X[:, 1:] - X[:, 1:].mean(0)) / X[:, 1:].std(0)
        b, *_ = np.linalg.lstsq(X, y, rcond=None)
        return b[1], b[2]
    obs_tm, obs_seq = coef(g.beta.values, g.tm.values, g[SEQ].values)
    nodes = sorted(set(g.a) | set(g.b))
    idx = {v: i for i, v in enumerate(nodes)}
    ia, ib = g.a.map(idx).values, g.b.map(idx).values
    M = np.full((len(nodes), len(nodes)), np.nan)
    M[ia, ib] = M[ib, ia] = g.beta.values
    c_tm = c_seq = valid = tries = 0
    while valid < n and tries < 20 * n:         # 湊滿 n 次有效置換，分母用有效次數
        tries += 1
        p = RNG.permutation(len(nodes))
        v = M[p[ia], p[ib]]
        ok = ~np.isnan(v)
        if ok.sum() < 30:
            continue
        valid += 1
        t, s_ = coef(v[ok], g.tm.values[ok], g[SEQ].values[ok])
        c_tm += abs(t) >= abs(obs_tm)
        c_seq += abs(s_) >= abs(obs_seq)
    assert valid == n
    return obs_tm, (c_tm + 1) / (valid + 1), obs_seq, (c_seq + 1) / (valid + 1)


def groups(g):
    """同源群＝以「共享註解」連起來的連通塊（cupin 7S＋11S 會在同一塊）。"""
    nodes = sorted(set(g.a) | set(g.b))
    parent = {v: v for v in nodes}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    for a, b in zip(g.a, g.b):
        parent[find(a)] = find(b)
    return {v: find(v) for v in nodes}


def main():
    g = homologous_core()
    rows = []
    r_tm = stats.spearmanr(g.tm, g.beta).statistic
    r_sq = stats.spearmanr(g[SEQ], g.beta).statistic
    r_xx = stats.spearmanr(g.tm, g[SEQ]).statistic
    b_tm, p_tm, b_sq, p_sq = incremental(g)
    print(f"同源對 n={len(g)}  ρ(TM)={r_tm:.3f}  ρ(seqid)={r_sq:.3f}  ρ(TM,seqid)={r_xx:.3f}")
    print(f"聯合模型（標準化秩）：TM β={b_tm:.2f} p={p_tm:.4f}；seqid β={b_sq:.2f} p={p_sq:.4f}")
    rows.append(dict(analysis="joint model", term="TM (rank, SD)", estimate=b_tm, p_perm=p_tm, n=len(g)))
    rows.append(dict(analysis="joint model", term="sequence identity (rank, SD)", estimate=b_sq, p_perm=p_sq, n=len(g)))
    rows.append(dict(analysis="collinearity", term="Spearman TM vs identity", estimate=r_xx, n=len(g)))
    rows.append(dict(analysis="full-length identity", term="rho identity",
                     estimate=stats.spearmanr(g.seqid, g.beta).statistic, n=len(g)))

    # 留一同源群
    comp = groups(g)
    g["grp"] = g.a.map(comp)
    sizes = g.grp.value_counts()
    loo = []
    for grp, k in sizes.items():
        h = g[g.grp != grp]
        loo.append(dict(left_out=grp, pairs_removed=int(k),
                        rho_tm=stats.spearmanr(h.tm, h.beta).statistic,
                        rho_seq=stats.spearmanr(h[SEQ], h.beta).statistic))
    loo = pd.DataFrame(loo)
    loo.to_csv(OUT / "robust_leave_one_group.csv", index=False)
    print(f"同源群 {len(sizes)} 個；留一後 ρ(TM) 範圍 {loo.rho_tm.min():.3f}–{loo.rho_tm.max():.3f}，"
          f"ρ(seqid) {loo.rho_seq.min():.3f}–{loo.rho_seq.max():.3f}")
    print(loo.sort_values("pairs_removed", ascending=False).head(6).to_string(index=False))
    rows.append(dict(analysis="leave one homology group out", term="rho TM min", estimate=loo.rho_tm.min(), n=len(sizes)))
    rows.append(dict(analysis="leave one homology group out", term="rho TM max", estimate=loo.rho_tm.max(), n=len(sizes)))

    # 分晶片
    bc = pd.read_csv(OUT / "cosens_by_chip.csv")
    bc = bc[bc.method == "mle"]
    for chip, c in bc.groupby("chip"):
        m = g.merge(c[["a", "b", "beta"]], on=["a", "b"], suffixes=("", "_chip"))
        r = stats.spearmanr(m.tm, m.beta_chip).statistic
        r2 = stats.spearmanr(m[SEQ], m.beta_chip).statistic
        print(f"{chip}: n={len(m)} ρ(TM)={r:.3f} ρ(seqid)={r2:.3f}")
        rows.append(dict(analysis="by chip", term=f"{chip} rho TM", estimate=r, n=len(m)))
        rows.append(dict(analysis="by chip", term=f"{chip} rho identity", estimate=r2, n=len(m)))

    # 門檻表
    bins = [0, 0.4, 0.5, 0.6, 0.7, 0.8, 1.01]
    g["tm_bin"] = pd.cut(g.tm, bins, right=False)
    t = g.groupby("tm_bin", observed=True).agg(
        n=("beta", "size"), median_logOR=("beta", "median"),
        q1=("beta", lambda x: x.quantile(.25)), q3=("beta", lambda x: x.quantile(.75)),
        pct_z_gt3=("z", lambda x: 100 * (x > 3).mean())).reset_index()
    t["tm_bin"] = t.tm_bin.astype(str)
    t.to_csv(OUT / "tm_bins_homologous.csv", index=False)
    print(t.to_string(index=False))
    # 共同暴露：同一目（order）的來源物種常一起吃／一起吸入（堅果、禾本科花粉），
    # 限制在來源物種屬於不同目的對，看結構的分級關係是否仍在
    lin = lineage_order()
    g["order_a"], g["order_b"] = g.a.map(lin), g.b.map(lin)
    h = g[g.order_a.notna() & g.order_b.notna() & (g.order_a != g.order_b)]
    r = stats.spearmanr(h.tm, h.beta).statistic
    r2 = stats.spearmanr(h[SEQ], h.beta).statistic
    print(f"來源物種不同目: n={len(h)} ρ(TM)={r:.3f} ρ(seqid)={r2:.3f}")
    rows.append(dict(analysis="different taxonomic order", term="rho TM", estimate=r, n=len(h)))
    rows.append(dict(analysis="different taxonomic order", term="rho identity", estimate=r2, n=len(h)))
    pd.DataFrame(rows).to_csv(OUT / "robust_summary.csv", index=False)


def lineage_order():
    """成分 → 來源物種所屬的目（UniProt lineage 中以 -ales/-ptera/-formes 等結尾者取第一個）。"""
    import json
    m = pd.read_csv(OUT / "component_uniprot.csv")
    m = m[(m.status == "ok") & m.accession.notna()].drop_duplicates("component")
    out = {}
    suf = ("ales", "ptera", "formes", "ida", "poda", "ina")
    for comp, acc in zip(m.component, m.accession):
        f = homology.CACHE / f"{acc}.json"
        if not f.exists():
            continue
        lin = json.loads(f.read_text()).get("organism", {}).get("lineage", [])
        o = [x for x in lin if (x.endswith(suf) or x in MAMMAL_ORDERS)
             and x not in ("Acari", "Arachnida")]
        out[comp] = o[-1] if o else None
    return out


if __name__ == "__main__":
    main()
