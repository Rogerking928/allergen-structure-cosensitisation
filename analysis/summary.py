"""稿件用的所有數字集中在這裡算一次，寫到 out/summary.json；正文只引用這個檔。"""
import json
import numpy as np
import pandas as pd
from scipy import stats
from common import load, ROOT
import structure_cosens as S

OUT = ROOT / "out"


def main():
    d, meta, comps = load()
    s = {}
    s["n_patients"] = len(d)
    s["chips"] = d.Chip_Type.value_counts().to_dict()
    s["n_component_columns"] = len(comps)
    m = pd.read_csv(OUT / "component_uniprot.csv")
    s["n_extract_or_ccd"] = int((m.status == "not_single_component").sum())
    s["n_no_accession"] = int((m.status == "no_accession").sum())
    s["n_mapped"] = int(m[m.status == "ok"].component.nunique())

    bc = pd.read_csv(OUT / "cosens_by_chip.csv")
    s["pairs_estimated_by_chip"] = bc.groupby("chip").size().to_dict()
    pooled = pd.read_csv(OUT / "cosens_pooled.csv")
    s["pairs_pooled"] = len(pooled)
    cru = pd.read_csv(OUT / "cosens_crude_pooled.csv")
    mm = pooled.merge(cru, on=["a", "b"])
    s["median_crude_logOR"] = float(mm.beta_crude.median())
    s["median_adjusted_logOR"] = float(mm.beta.median())

    p = pd.read_csv(OUT / "pairs_structure_cosens.csv")
    s["pairs_with_structure"] = len(p)
    s["excl_same_species"] = int(p.same_species.sum())
    q = p[~p.same_species]
    s["excl_subchain"] = int(q.subchain.sum()); q = q[~q.subchain]
    short = (q.len_a < S.MIN_LEN) | (q.len_b < S.MIN_LEN)
    s["excl_short"] = int(short.sum()); q = q[~short]
    s["excl_unannotated"] = int((~q.annotated).sum()); q = q[q.annotated]
    s["pairs_core"] = len(q)
    hom, non = q[q.homologous], q[~q.homologous]
    s["pairs_homologous"] = len(hom)
    s["pairs_nonhomologous"] = len(non)
    s["components_in_core"] = len(set(q.a) | set(q.b))
    s["nonhomologous_tm_ge_05"] = int((non.tm >= S.TM_HI).sum())
    s["nonhomologous_tm_max"] = float(non.tm.max())
    s["nonhomologous_tm_p99"] = float(non.tm.quantile(.99))
    ss = p[p.same_species & ~p.subchain]
    for name, g in [("same_species", ss), ("homologous", hom), ("nonhomologous", non)]:
        s[f"{name}_n"] = len(g)
        s[f"{name}_median"] = float(g.beta.median())
        s[f"{name}_q1"] = float(g.beta.quantile(.25))
        s[f"{name}_q3"] = float(g.beta.quantile(.75))
        s[f"{name}_pct_z3"] = float(100 * (g.z > 3).mean())

    lo = hom[hom.tm < S.TM_HI]
    fams = {"2S albumin (prolamin)", "nsLTP"}
    s["hom_tm_lt05_n"] = len(lo)
    s["hom_tm_lt05_2S_vs_LTP"] = int(((lo.fam_a != lo.fam_b) & lo.fam_a.isin(fams) & lo.fam_b.isin(fams)).sum())
    s["hom_tm_lt05_codex"] = int((lo.seqid80 >= S.SEQ_HI).sum())
    hi = hom[hom.tm >= S.TM_HI]
    s["hom_tm_ge05_n"] = len(hi)
    s["hom_tm_ge05_codex_neg"] = int((hi.seqid80 < S.SEQ_HI).sum())
    s["hom_codex_neg_n"] = int((hom.seqid80 < S.SEQ_HI).sum())
    s["hom_codex_neg_median_ge05"] = float(hom[(hom.seqid80 < S.SEQ_HI) & (hom.tm >= S.TM_HI)].beta.median())
    s["hom_codex_neg_median_lt05"] = float(hom[(hom.seqid80 < S.SEQ_HI) & (hom.tm < S.TM_HI)].beta.median())
    s["hom_codex_neg_ge05_n"] = int(((hom.seqid80 < S.SEQ_HI) & (hom.tm >= S.TM_HI)).sum())
    ms = pd.read_csv(OUT / "structure_mantel.csv").set_index("stratum")
    s["mantel"] = ms.reset_index().to_dict(orient="records")
    s["robust"] = pd.read_csv(OUT / "robust_summary.csv").to_dict(orient="records")
    s["tm_bins"] = pd.read_csv(OUT / "tm_bins_homologous.csv").to_dict(orient="records")
    loo = pd.read_csv(OUT / "robust_leave_one_group.csv")
    s["n_homology_groups"] = len(loo)
    json.dump(s, open(OUT / "summary.json", "w"), indent=1, default=float)
    print(json.dumps(s, indent=1, default=float)[:6000])


if __name__ == "__main__":
    main()
