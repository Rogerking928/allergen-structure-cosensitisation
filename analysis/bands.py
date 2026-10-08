"""TM-score 分帶的勝算比中位數（含節點自助法信賴區間）與排列檢定的虛無分布。

Table 2 要報的是 median OR 與 **95% CI**（IQR 描述離散度，不是推論），
趨勢 P 用的是同源對內 Spearman 的排列檢定 P（與 structure_cosens 同一個做法）。
虛無分布另外存起來，供補充圖畫出來——審稿人看到排列檢定時會想知道虛無長什麼樣。
"""
import json
import numpy as np
import pandas as pd
from scipy import stats
from common import ROOT
import structure_cosens as S
import discrim

OUT = ROOT / "out"
RNG = np.random.default_rng(20261007)
BINS = [0, 0.4, 0.5, 0.6, 0.7, 0.8, 1.01]
LABEL = {"[0.0, 0.4)": "< 0.40", "[0.4, 0.5)": "0.40 to < 0.50", "[0.5, 0.6)": "0.50 to < 0.60",
         "[0.6, 0.7)": "0.60 to < 0.70", "[0.7, 0.8)": "0.70 to < 0.80", "[0.8, 1.01)": "≥ 0.80"}
N_BOOT = 2000


def median_ci(g, rng, n=N_BOOT):
    vals = []
    for _ in range(n):
        w = discrim.node_boot_weights(g, rng)
        if w.sum() < 5:
            continue
        idx = np.repeat(np.arange(len(g)), w.astype(int))
        vals.append(np.median(g.beta.values[idx]))
    return np.percentile(vals, [2.5, 97.5])


def nulls(g, x="tm", n=10000):
    """排列檢定的虛無分布（打亂過敏原標籤，不是打亂對）。"""
    nodes = sorted(set(g.a) | set(g.b))
    idx = {v: i for i, v in enumerate(nodes)}
    ia, ib = g.a.map(idx).values, g.b.map(idx).values
    M = np.full((len(nodes), len(nodes)), np.nan)
    M[ia, ib] = M[ib, ia] = g.beta.values
    rng = np.random.default_rng(20261003)
    out = []
    tries = 0
    while len(out) < n and tries < 20 * n:      # 湊滿 n 次有效置換（見 structure_cosens.mantel）
        tries += 1
        p = rng.permutation(len(nodes))
        v = M[p[ia], p[ib]]
        ok = ~np.isnan(v)
        if ok.sum() < 30:
            continue
        out.append(stats.spearmanr(g[x].values[ok], v[ok]).statistic)
    assert len(out) == n
    return np.array(out)


def main():
    c = discrim.core_pairs()
    hom, non = c[c.homologous], c[~c.homologous]
    hom = hom.copy()
    hom["bin"] = pd.cut(hom.tm, BINS, right=False).astype(str).map(LABEL)
    rng = np.random.default_rng(RNG.integers(1 << 32))
    rows = []
    for lab in LABEL.values():
        g = hom[hom["bin"] == lab]
        lo, hi = median_ci(g, rng)
        ex = g.assign(d=(g.beta - g.beta.median()).abs()).sort_values("d").iloc[0]
        rows.append(dict(band=lab, n=len(g), median_or=np.exp(g.beta.median()),
                         lo=np.exp(lo), hi=np.exp(hi), pct_z3=100 * (g.z > 3).mean(),
                         example=f"{ex.a.replace('_RUO', '').replace('_', ' ')} – "
                                 f"{ex.b.replace('_RUO', '').replace('_', ' ')}",
                         example_or=np.exp(ex.beta)))
    t = pd.DataFrame(rows)
    t.to_csv(OUT / "tm_bands_ci.csv", index=False)
    print(t.to_string(index=False))

    nh = nulls(hom)
    nn = nulls(non)
    np.save(OUT / "null_homologous.npy", nh)
    np.save(OUT / "null_nonhomologous.npy", nn)
    obs_h = stats.spearmanr(hom.tm, hom.beta).statistic
    obs_n = stats.spearmanr(non.tm, non.beta).statistic
    res = dict(obs_homologous=float(obs_h), obs_nonhomologous=float(obs_n),
               null_homologous_p975=float(np.percentile(nh, 97.5)),
               null_nonhomologous_p975=float(np.percentile(nn, 97.5)),
               trend_p_homologous=float((np.sum(np.abs(nh) >= abs(obs_h)) + 1) / (len(nh) + 1)))
    (OUT / "bands.json").write_text(json.dumps(res, indent=1))
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
