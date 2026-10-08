"""2026-10-08 第三方核對報告要求的補強分析。

逐項對應審閱意見：
1. 多重比較：18,848 對的 BH q 值（寫進 Supplementary Data S1）。
2. AUC 比較不可只看信賴區間重疊：DeLong 配對檢定（同一組對），
   另加**節點自助法的 AUC 差值**——DeLong 假設觀察獨立，而成分對不獨立，兩個都報。
3. z > 3 是任意的：z > 2 與 z > 4 重算。
4. 長鏈正規化的代價：改用短鏈正規化（tm_max）重算 ρ 與「只有結構達標」那一格。
5. 平台異質性：每一對的 Q 與 I²。
6. 二分化損失資訊：改用連續 IgE 的**偏 Spearman**（對廣度取偏相關）重算。
7. 廣度校正的 collider 疑慮：廣度改用「與兩端都不同家族」的成分計算，重估全部核心對。
8. 核對 Discussion 的 Arg r 1 人數（32/1,139、其中 8 人同時對 Can f 1 陽性）。
"""
import itertools
import json

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy import stats
from sklearn.metrics import roc_auc_score

import cosens
import discrim
import families
import structure_cosens as S
from common import ROOT, load

OUT = ROOT / "out"
RNG = np.random.default_rng(20261008)


# ----------------------------------------------------------------- 1. BH q 值
def qvalues():
    p = pd.read_csv(OUT / "cosens_pooled.csv")
    p["p_value"] = 2 * stats.norm.sf(np.abs(p.z))
    o = np.argsort(p.p_value.values)
    ranked = p.p_value.values[o]
    q = ranked * len(p) / (np.arange(len(p)) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty(len(p))
    out[o] = np.minimum(q, 1.0)
    p["q_value"] = out
    p[["a", "b", "p_value", "q_value"]].to_csv(OUT / "qvalues.csv", index=False)
    c = discrim.core_pairs().merge(p[["a", "b", "q_value"]], on=["a", "b"])
    sub = c[c.struct & ~c.codex]
    return dict(n=len(p), q05=int((p.q_value < 0.05).sum()),
                core_q05=int((c.q_value < 0.05).sum()),
                strong_q05=int(((c.z > 3) & (c.q_value < 0.05)).sum()),
                struct_only_n=len(sub), struct_only_q05=int((sub.q_value < 0.05).sum()),
                struct_only_strong=int((sub.z > 3).sum()))


# --------------------------------------------- 2. DeLong ＋ 節點自助法的 AUC 差
def _midrank(x):
    j = np.argsort(x)
    z = x[j]
    n = len(x)
    t = np.zeros(n)
    i = 0
    while i < n:
        k = i
        while k < n and z[k] == z[i]:
            k += 1
        t[i:k] = 0.5 * (i + k - 1) + 1
        i = k
    out = np.empty(n)
    out[j] = t
    return out


def delong(y, xs):
    """兩條以上相關 ROC 曲線的 DeLong 共變異數（Sun & Xu 2014 的快速版）。"""
    y = np.asarray(y).astype(bool)
    pos = np.array([np.asarray(x)[y] for x in xs])
    neg = np.array([np.asarray(x)[~y] for x in xs])
    m, n = pos.shape[1], neg.shape[1]
    k = len(xs)
    tz = np.empty((k, m)), np.empty((k, n))
    v01, v10 = np.empty((k, m)), np.empty((k, n))
    aucs = np.empty(k)
    for i in range(k):
        tx, ty = _midrank(pos[i]), _midrank(neg[i])
        txy = _midrank(np.concatenate([pos[i], neg[i]]))
        aucs[i] = (txy[:m].sum() - m * (m + 1) / 2) / (m * n)
        v01[i] = (txy[:m] - tx) / n
        v10[i] = 1 - (txy[m:] - ty) / m
    s = np.cov(v01) / m + np.cov(v10) / n
    s = np.atleast_2d(s)
    d = aucs[0] - aucs[1]
    var = s[0, 0] + s[1, 1] - 2 * s[0, 1]
    z = d / np.sqrt(var) if var > 0 else np.nan
    return dict(auc1=aucs[0], auc2=aucs[1], diff=d, se=np.sqrt(var),
                z=z, p=2 * stats.norm.sf(abs(z)))


def auc_diff_boot(g, c1="tm", c2="seqid80", n=2000):
    """AUC 差值的節點自助法區間（成分抽後放回，一對的權重＝兩端次數相乘）。"""
    y, x1, x2 = g.strong.values, g[c1].values, g[c2].values
    rng = np.random.default_rng(RNG.integers(1 << 32))
    vals = []
    for _ in range(n):
        w = discrim.node_boot_weights(g, rng)
        if w.sum() == 0 or len(np.unique(y[w > 0])) < 2:
            continue
        vals.append(roc_auc_score(y, x1, sample_weight=w)
                    - roc_auc_score(y, x2, sample_weight=w))
    lo, hi = np.percentile(vals, [2.5, 97.5])
    return dict(diff=roc_auc_score(y, x1) - roc_auc_score(y, x2), lo=lo, hi=hi, n_boot=len(vals))


# ------------------------------------------------------- 3. z 門檻的敏感度分析
def z_thresholds(c):
    rows = []
    for z in (2.0, 3.0, 4.0):
        s = c.z > z
        cells = {}
        for st in (True, False):
            for cx in (True, False):
                g = c[(c.struct == st) & (c.codex == cx)]
                cells[f"{'S' if st else 's'}{'C' if cx else 'c'}"] = \
                    dict(n=len(g), strong=int((g.z > z).sum()))
        rows.append(dict(z=z, n_strong=int(s.sum()), pct=100 * s.mean(),
                         auc_tm=roc_auc_score(s, c.tm), auc_seq=roc_auc_score(s, c.seqid80),
                         struct_only_strong=cells["Sc"]["strong"], struct_only_n=cells["Sc"]["n"],
                         neither_pct=100 * cells["sc"]["strong"] / cells["sc"]["n"]))
    return rows


# ------------------------------------------------- 4. 短鏈正規化的敏感度分析
def short_chain(c):
    hom = c[c.homologous]
    out = dict(
        rho_long=stats.spearmanr(hom.tm, hom.beta).statistic,
        rho_short=stats.spearmanr(hom.tm_max, hom.beta).statistic,
        rho_sym=stats.spearmanr(hom.tm_sym, hom.beta).statistic,
        nonhom_max_long=float(c[~c.homologous].tm.max()),
        nonhom_max_short=float(c[~c.homologous].tm_max.max()),
        nonhom_ge05_short=int((c[~c.homologous].tm_max >= S.TM_HI).sum()),
    )
    for name, col in (("long", "tm"), ("short", "tm_max")):
        st = c[col] >= S.TM_HI
        only = c[st & ~c.codex]
        out[f"struct_only_{name}"] = len(only)
        out[f"struct_only_{name}_strong"] = int((only.z > 3).sum())
        out[f"seq_only_{name}"] = int((~st & c.codex).sum())
        out[f"auc_{name}"] = roc_auc_score(c.z > 3, c[col])
    return out


# ------------------------------------------------------------ 5. 平台異質性
def heterogeneity():
    bc = pd.read_csv(OUT / "cosens_by_chip.csv")
    bc = bc[(bc.method == "mle") & bc.se.notna()].copy()
    bc["w"] = 1 / bc.se ** 2
    rows = []
    for (a, b), g in bc.groupby(["a", "b"]):
        if len(g) < 2:
            continue
        mu = np.average(g.beta, weights=g.w)
        q = float((g.w * (g.beta - mu) ** 2).sum())
        df = len(g) - 1
        rows.append(dict(a=a, b=b, k=len(g), Q=q, df=df,
                         I2=max(0.0, 100 * (q - df) / q) if q > 0 else 0.0,
                         p_het=float(stats.chi2.sf(q, df))))
    h = pd.DataFrame(rows)
    h.to_csv(OUT / "heterogeneity.csv", index=False)
    core = discrim.core_pairs()[["a", "b"]].merge(h, on=["a", "b"])
    return dict(n_pairs_multi=len(h), median_I2=float(h.I2.median()),
                pct_I2_gt50=float(100 * (h.I2 > 50).mean()),
                pct_het_p05=float(100 * (h.p_het < 0.05).mean()),
                core_n=len(core), core_median_I2=float(core.I2.median()),
                core_pct_I2_gt50=float(100 * (core.I2 > 50).mean()))


# ------------------------------------- 6. 連續 IgE（對廣度取偏相關的 Spearman）
def _partial_spearman(x, y, z):
    r = np.column_stack([stats.rankdata(x), stats.rankdata(y), stats.rankdata(z)])
    r = (r - r.mean(0)) / r.std(0)
    zz = np.column_stack([np.ones(len(r)), r[:, 2]])
    bx = np.linalg.lstsq(zz, r[:, 0], rcond=None)[0]
    by = np.linalg.lstsq(zz, r[:, 1], rcond=None)[0]
    ex, ey = r[:, 0] - zz @ bx, r[:, 1] - zz @ by
    return float(np.corrcoef(ex, ey)[0, 1])


def continuous(core):
    d, meta, comps = load()
    need = set(core.a) | set(core.b)
    vals = {}
    for chip in d.Chip_Type.unique():
        y, keep = cosens.strat_table(d, comps, chip)
        g = d[d.Chip_Type == chip]
        x = g[comps].replace(-1, np.nan)
        x = x.loc[y.index, [c for c in keep]]
        breadth_all = (x >= cosens.POS).sum(axis=1).values
        cols = {c: x[c].values for c in keep if c in need}
        vals[chip] = (cols, breadth_all)
    rows = []
    for t in core.itertuples():
        rs, ws = [], []
        for chip, (cols, br) in vals.items():
            if t.a not in cols or t.b not in cols:
                continue
            xa, xb = cols[t.a], cols[t.b]
            m = ~(np.isnan(xa) | np.isnan(xb))
            if m.sum() < 50:
                continue
            # 廣度扣掉這一對自己，避免把兩端算進條件變數
            br_pair = br[m] - (xa[m] >= cosens.POS) - (xb[m] >= cosens.POS)
            if len(np.unique(br_pair)) < 3:
                continue
            r = _partial_spearman(xa[m], xb[m], br_pair)
            if np.isfinite(r):
                rs.append(np.arctanh(np.clip(r, -0.999, 0.999)))
                ws.append(m.sum() - 4)
        if rs:
            rows.append(dict(a=t.a, b=t.b, rho_cont=np.tanh(np.average(rs, weights=ws)),
                             k=len(rs)))
    r = pd.DataFrame(rows)
    r.to_csv(OUT / "continuous_partial_spearman.csv", index=False)
    m = core.merge(r, on=["a", "b"])
    hom, non = m[m.homologous], m[~m.homologous]
    return dict(n=len(m),
                rho_tm_vs_continuous_homologous=stats.spearmanr(hom.tm, hom.rho_cont).statistic,
                rho_seq_vs_continuous_homologous=stats.spearmanr(hom.seqid80,
                                                                 hom.rho_cont).statistic,
                rho_tm_vs_continuous_nonhomologous=stats.spearmanr(non.tm, non.rho_cont).statistic,
                agreement_with_binary=stats.spearmanr(m.beta, m.rho_cont).statistic)


# --------------------------- 7. 廣度改用「與兩端都不同家族」的成分（collider）
def disjoint_breadth(core):
    d, meta, comps = load()
    fam = families.table()
    fam = fam[fam.status == "ok"].drop_duplicates("component").set_index("component").family
    need = set(core.a) | set(core.b)
    per_chip = {}
    for chip in d.Chip_Type.unique():
        y, keep = cosens.strat_table(d, comps, chip)
        arr = y[keep].values
        idx = {c: i for i, c in enumerate(keep)}
        famv = np.array([str(fam.get(c, "NA")) for c in keep])
        per_chip[chip] = (arr, idx, famv)
    rows = []
    for t in core.itertuples():
        fa, fb = str(fam.get(t.a, "NA")), str(fam.get(t.b, "NA"))
        for chip, (arr, idx, famv) in per_chip.items():
            if t.a not in idx or t.b not in idx:
                continue
            ya, yb = arr[:, idx[t.a]], arr[:, idx[t.b]]
            m = ~(np.isnan(ya) | np.isnan(yb))
            if m.sum() < 50 or ((ya[m] == 1) & (yb[m] == 1)).sum() < cosens.MIN_CELL:
                continue
            out_fam = (famv != fa) & (famv != fb)      # 與兩端都不同家族的成分
            br = np.nansum(arr[np.ix_(m, np.where(out_fam)[0])], axis=1)
            beta, se, how = cosens.fit_pair(ya[m], yb[m], br)
            if how == "mle" and np.isfinite(se):
                rows.append(dict(a=t.a, b=t.b, chip=chip, beta=beta, se=se))
    r = pd.DataFrame(rows)
    r["w"] = 1 / r.se ** 2
    g = r.groupby(["a", "b"])
    pooled = pd.DataFrame({"beta_dj": g.apply(lambda x: np.average(x.beta, weights=x.w),
                                              include_groups=False),
                           "se_dj": 1 / np.sqrt(g.w.sum())}).reset_index()
    pooled.to_csv(OUT / "cosens_disjoint_breadth.csv", index=False)
    m = core.merge(pooled, on=["a", "b"])
    hom, non = m[m.homologous], m[~m.homologous]
    return dict(n=len(m), median_beta_all=float(m.beta_dj.median()),
                median_beta_nonhom=float(non.beta_dj.median()),
                median_or_nonhom=float(np.exp(non.beta_dj.median())),
                rho_tm_homologous=stats.spearmanr(hom.tm, hom.beta_dj).statistic,
                rho_seq_homologous=stats.spearmanr(hom.seqid80, hom.beta_dj).statistic,
                rho_tm_nonhomologous=stats.spearmanr(non.tm, non.beta_dj).statistic,
                agreement_with_main=stats.spearmanr(m.beta, m.beta_dj).statistic)


# ------------------------------------------------- 8. 核對 Arg r 1 的逐人數字
def arg_r_1():
    d, meta, comps = load()
    out = {}
    for chip in d.Chip_Type.unique():
        g = d[d.Chip_Type == chip]
        if "Arg_r_1" not in g or g["Arg_r_1"].replace(-1, np.nan).notna().mean() < 0.8:
            continue
        x = g[["Arg_r_1", "Can_f_1"]].replace(-1, np.nan)
        pa = x.Arg_r_1 >= cosens.POS
        pb = x.Can_f_1 >= cosens.POS
        out[chip] = dict(patients=int(len(g)), arg_r_1_positive=int(pa.sum()),
                         both_positive=int((pa & pb).sum()),
                         can_f_1_positive=int(pb.sum()))
    return out


def main():
    core = discrim.core_pairs()
    res = {}
    print("1 q 值", flush=True)
    res["fdr"] = qvalues()
    print("2 AUC 比較", flush=True)
    res["delong_all"] = delong(core.strong.values, [core.tm.values, core.seqid80.values])
    hom = core[core.homologous]
    res["delong_homologous"] = delong(hom.strong.values, [hom.tm.values, hom.seqid80.values])
    res["auc_diff_boot_all"] = auc_diff_boot(core)
    res["auc_diff_boot_homologous"] = auc_diff_boot(hom)
    print("3 z 門檻", flush=True)
    res["z_thresholds"] = z_thresholds(core)
    print("4 短鏈正規化", flush=True)
    res["normalisation"] = short_chain(core)
    print("5 異質性", flush=True)
    res["heterogeneity"] = heterogeneity()
    print("8 Arg r 1", flush=True)
    res["arg_r_1"] = arg_r_1()
    print("6 連續 IgE（慢）", flush=True)
    res["continuous"] = continuous(core)
    print("7 不相交廣度（最慢）", flush=True)
    res["disjoint_breadth"] = disjoint_breadth(core)
    (OUT / "extra.json").write_text(json.dumps(res, indent=1, default=float))
    print(json.dumps(res, indent=1, default=float))


if __name__ == "__main__":
    main()
