"""AlphaFold 模型信心（pLDDT）會不會汙染結構相似度。

審稿人一定會問：低信心的模型（無序區、片段）算出來的 TM-score 可不可信。
AlphaFold 把 pLDDT 放在 PDB 的 B-factor 欄；本檔逐個 accession 算 CA 原子的平均 pLDDT，
再做一個「兩端模型都 ≥ 70」的敏感度分析。

70 是 AlphaFold 官方的分界（> 70 視為骨架大致可信）。
"""
import json
import numpy as np
import pandas as pd
from scipy import stats
from common import ROOT
import families
import structure_cosens as S

OUT = ROOT / "out"
AFDB = ROOT / "data" / "afdb"
CUT = 70.0


def mean_plddt(acc):
    p = AFDB / f"{acc}.pdb"
    if not p.exists():
        return np.nan
    vals = []
    for line in p.read_text().splitlines():
        if line.startswith("ATOM") and line[12:16].strip() == "CA":
            vals.append(float(line[60:66]))
    return float(np.mean(vals)) if vals else np.nan


def by_component():
    m = families.table()
    m = m[(m.status == "ok") & m.accession.notna()]
    rows = []
    for comp, g in m.groupby("component"):
        v = [mean_plddt(a) for a in g.accession]
        v = [x for x in v if not np.isnan(x)]
        if v:
            rows.append(dict(component=comp, plddt=float(np.mean(v))))
    return pd.DataFrame(rows)


def main():
    cp = by_component()
    cp.to_csv(OUT / "component_plddt.csv", index=False)
    lut = cp.set_index("component").plddt

    d = pd.read_csv(OUT / "pairs_structure_cosens.csv")
    core = d[(~d.same_species) & (~d.subchain) & (d.len_a >= S.MIN_LEN) & (d.len_b >= S.MIN_LEN)
             & d.annotated].copy()
    core["plddt_min"] = np.minimum(core.a.map(lut), core.b.map(lut))
    hom = core[core.homologous].dropna(subset=["plddt_min"])
    keep = hom[hom.plddt_min >= CUT]

    res = dict(cut=CUT, n_components=len(cp),
               median_plddt=float(cp.plddt.median()),
               q1=float(cp.plddt.quantile(.25)), q3=float(cp.plddt.quantile(.75)),
               n_components_below=int((cp.plddt < CUT).sum()),
               n_homologous=len(hom), n_homologous_kept=len(keep))
    r_all, p_all = S.mantel(hom)
    r_hi, p_hi = S.mantel(keep)
    res["rho_all"], res["p_all"] = r_all, p_all
    res["rho_high_confidence"], res["p_high_confidence"] = r_hi, p_hi
    # pLDDT 本身與 TM-score 的關係（低信心模型是否系統性地算出較低的 TM）
    res["rho_plddt_tm"] = float(stats.spearmanr(hom.plddt_min, hom.tm).statistic)
    (OUT / "plddt.json").write_text(json.dumps(res, indent=1))
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
