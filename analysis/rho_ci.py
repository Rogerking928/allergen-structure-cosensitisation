"""敏感度分析那些 Spearman ρ 的信賴區間（節點自助法）。

原本只報點估計，補充圖看起來像沒做推論。重抽的是**過敏原**不是對，
理由與 discrim.py 相同：同一個成分出現在很多對裡。
"""
import json
import numpy as np
import pandas as pd
from scipy import stats
from common import ROOT
import discrim
import structure_cosens as S

OUT = ROOT / "out"
N_BOOT = 2000
RNG = np.random.default_rng(20261007)


def rho_ci(g, x="tm", y="beta", n=N_BOOT):
    obs = stats.spearmanr(g[x], g[y]).statistic
    rng = np.random.default_rng(RNG.integers(1 << 32))
    vals = []
    for _ in range(n):
        w = discrim.node_boot_weights(g, rng)
        if (w > 0).sum() < 20:
            continue
        idx = np.repeat(np.arange(len(g)), w.astype(int))
        if len(idx) < 20:
            continue
        vals.append(stats.spearmanr(g[x].values[idx], g[y].values[idx]).statistic)
    lo, hi = np.percentile(vals, [2.5, 97.5])
    return dict(n=len(g), rho=float(obs), lo=float(lo), hi=float(hi))


def main():
    import robust
    c = discrim.core_pairs()
    hom, non = c[c.homologous], c[~c.homologous]
    cp = pd.read_csv(OUT / "component_plddt.csv").set_index("component").plddt
    bc = pd.read_csv(OUT / "cosens_by_chip.csv")
    bc = bc[bc.method == "mle"]
    lin = robust.lineage_order()

    rows = {}
    rows["homologous"] = rho_ci(hom)
    for chip, lab in [("ISAC_V1", "isac_v1"), ("ISAC_V2", "isac_v2"), ("ALEX", "alex")]:
        m = hom.drop(columns=["beta"]).merge(
            bc[bc.chip == chip][["a", "b", "beta"]], on=["a", "b"])
        rows[lab] = rho_ci(m)
    h = hom.assign(oa=hom.a.map(lin), ob=hom.b.map(lin))
    rows["different_order"] = rho_ci(h[h.oa.notna() & h.ob.notna() & (h.oa != h.ob)])
    pl = hom.assign(p=np.minimum(hom.a.map(cp), hom.b.map(cp)))
    rows["plddt70"] = rho_ci(pl[pl.p >= 70])
    rows["below_codex"] = rho_ci(hom[hom.seqid80 < S.SEQ_HI])
    rows["nonhomologous"] = rho_ci(non)
    rows["nonhomologous_nofrag"] = rho_ci(non[~non.frag])
    rows["nonhomologous_curated"] = rho_ci(c[~c.same_family])
    (OUT / "rho_ci.json").write_text(json.dumps(rows, indent=1))
    for k, v in rows.items():
        print(f"{k:24s} n={v['n']:5d}  rho {v['rho']:.3f} ({v['lo']:.3f} to {v['hi']:.3f})")


if __name__ == "__main__":
    main()
