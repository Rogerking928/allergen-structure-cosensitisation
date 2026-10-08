"""結構相似度能不能預測共同致敏——本專案的主結果。

分析單位是**成分對**。對不是獨立的（同一個成分出現在很多對裡，而且家族成群），
所以顯著性一律用 Mantel 式的排列檢定：打亂的是**過敏原的標籤**，不是對。

主分析排除：
- 同物種的對（共同暴露，不是交叉反應）
- 同一前驅物切出來的結構域（Hev b 6.01／6.02，結構相似度是人工的）
重點放在**跨家族**的對：序列上看不出關係、但結構像的那些，是否也常一起陽性。
"""
import numpy as np
import pandas as pd
from scipy import stats
from common import ROOT
import families
import homology

OUT = ROOT / "out"
TM_HI = 0.5          # 慣用的同摺疊門檻
MIN_LEN = 50         # 太短的片段模型不能比對：21 殘基的 Ole e 7 對什麼都 >0.9
SEQ_HI = 35.0        # FAO/WHO Codex：80 殘基窗內 identity ≥ 35% 視為可能交叉反應
RNG = np.random.default_rng(20261002)


def build():
    m = families.table()
    m = m[(m.status == "ok") & m.accession.notna()]
    # 多鏈（Fel d 1）：一個成分對到兩個 accession，相似度取鏈間最大值
    comp2acc = m.groupby("component").accession.apply(list).to_dict()
    meta = m.drop_duplicates("component").set_index("component")

    sim = pd.read_csv(OUT / "pair_similarity.csv")
    sim["tm_sym"] = sim.tm                      # foldseek alntmscore（敏感度分析用）
    # pair_seqid 的 seqid 是「局部比對對到的那一段」的 identity（短片段常是 100%），不是全長。
    # 全長 identity＝同一個局部比對的相同殘基數 ÷ 較長序列的長度（2026-10-08 審閱抓到）
    ps = pd.read_csv(OUT / "pair_seqid.csv")
    seqlen = m.drop_duplicates("accession").set_index("accession").sequence.str.len()
    ps["matches"] = ps.seqid * ps.alnlen / 100
    ps["seqid_full"] = 100 * ps.matches / np.maximum(ps.a.map(seqlen), ps.b.map(seqlen))
    full = {}
    for a, b, v in ps[["a", "b", "seqid_full"]].itertuples(index=False):
        full[(a, b)] = full[(b, a)] = v
    key = {}
    for a, b, tm, tmax, tmin, sid, s80 in sim[["a", "b", "tm", "tm_max", "tm_min",
                                               "seqid", "seqid_best80"]].itertuples(index=False):
        key[(a, b)] = key[(b, a)] = (tm, tmax, tmin, full.get((a, b), 0.0), s80, sid)

    p = pd.read_csv(OUT / "cosens_pooled.csv")
    rows = []
    for t in p.itertuples():
        A, B = comp2acc.get(t.a), comp2acc.get(t.b)
        if not A or not B:
            continue
        vals = [key[(x, y)] for x in A for y in B if (x, y) in key]
        if not vals:
            continue
        v = max(vals, key=lambda z: z[1])       # 鏈間取結構最像的一組
        ma, mb = meta.loc[t.a], meta.loc[t.b]
        rows.append(dict(a=t.a, b=t.b, beta=t.beta, se=t.se, z=t.z, k_chip=t.k_chip,
                         fam_a=ma.family, fam_b=mb.family, org_a=ma.organism, org_b=mb.organism,
                         frag=bool(ma.fragment or mb.fragment),
                         subchain=bool(str(ma.note) == "subchain" or str(mb.note) == "subchain"),
                         tm_sym=v[0], tm_max=v[1], tm=v[2], seqid=v[3], seqid80=v[4],
                         seqid_local=v[5],
                         len_a=ma.length, len_b=mb.length))
    d = pd.DataFrame(rows)
    d["same_species"] = d.org_a == d.org_b
    d["same_family"] = (d.fam_a == d.fam_b) & (d.fam_a != "Other / unclassified")
    ann = homology.table()
    d["annotated"] = d.a.map(lambda c: bool(ann.get(c))) & d.b.map(lambda c: bool(ann.get(c)))
    d["homologous"] = [homology.homologous(ann, a, b) for a, b in zip(d.a, d.b)]
    return d


def mantel(d, x="tm", y="beta", n=10000):
    """排列檢定：把過敏原標籤重排，重算同一組對的統計量。"""
    nodes = sorted(set(d.a) | set(d.b))
    obs = stats.spearmanr(d[x], d[y]).statistic
    idx = {v: i for i, v in enumerate(nodes)}
    ia, ib = d.a.map(idx).values, d.b.map(idx).values
    M = np.full((len(nodes), len(nodes)), np.nan)
    M[ia, ib] = M[ib, ia] = d[y].values
    # 重排後可比的對不足 30 的那次不算；一直抽到湊滿 n 次「有效」置換，P 的分母用有效次數
    # （同源對很稀疏，原本約三成被跳過、分母卻仍寫 n，P 會偏小）
    # 每次呼叫用同一個種子、與 bands.nulls 相同：同源／非同源的 P 才會和補充圖 S3、Source data 的虛無分布逐次一致
    rng = np.random.default_rng(20261003)
    cnt = valid = tries = 0
    while valid < n and tries < 20 * n:
        tries += 1
        perm = rng.permutation(len(nodes))
        v = M[perm[ia], perm[ib]]
        ok = ~np.isnan(v)
        if ok.sum() < 30:
            continue
        valid += 1
        r = stats.spearmanr(d[x].values[ok], v[ok]).statistic
        cnt += abs(r) >= abs(obs)
    assert valid == n, f"有效置換只有 {valid} 次"
    return obs, (cnt + 1) / (valid + 1)


def describe(name, g):
    if len(g) < 5:
        return dict(stratum=name, n=len(g))
    return dict(stratum=name, n=len(g),
                median_beta=g.beta.median(),
                q1=g.beta.quantile(.25), q3=g.beta.quantile(.75),
                pct_z_gt3=100 * (g.z > 3).mean())


def main():
    d = build()
    d.to_csv(OUT / "pairs_structure_cosens.csv", index=False)
    core = d[(~d.same_species) & (~d.subchain) & (d.len_a >= MIN_LEN) & (d.len_b >= MIN_LEN)
             & d.annotated]
    # 主分析的「跨家族」改用 Pfam／SUPFAM 同源（homology.py）；手寫家族標籤只留作敏感度分析
    hom = core[core.homologous]
    cross = core[~core.homologous]
    print(f"對數：全部 {len(d)}、排除同物種／同前驅物／<{MIN_LEN} aa 片段／無註解後 {len(core)}、"
          f"同源 {len(hom)}、非同源 {len(cross)}")

    out = []
    out.append(describe("All cross-species pairs", core))
    out.append(describe("Homologous (shared Pfam or SCOP superfamily)", hom))
    out.append(describe("Non-homologous", cross))
    out.append(describe(f"Non-homologous, TM ≥ {TM_HI}", cross[cross.tm >= TM_HI]))
    out.append(describe(f"Non-homologous, TM < {TM_HI}", cross[cross.tm < TM_HI]))
    out.append(describe(f"Homologous, TM ≥ {TM_HI}", hom[hom.tm >= TM_HI]))
    out.append(describe(f"Homologous, TM < {TM_HI}", hom[hom.tm < TM_HI]))
    out.append(describe(f"Homologous, 80-aa identity ≥ {SEQ_HI}%", hom[hom.seqid80 >= SEQ_HI]))
    out.append(describe(f"Homologous, 80-aa identity < {SEQ_HI}%", hom[hom.seqid80 < SEQ_HI]))
    s = pd.DataFrame(out)
    s.to_csv(OUT / "structure_strata.csv", index=False)
    print(s.to_string(index=False))

    res = []
    for name, g, x in [("all cross-species", core, "tm"),
                       ("non-homologous", cross, "tm"),
                       ("homologous", hom, "tm"),
                       ("homologous, sequence identity", hom, "seqid80"),
                       ("homologous, full-length identity", hom, "seqid"),
                       ("homologous, below Codex sequence criterion", hom[hom.seqid80 < SEQ_HI], "tm"),
                       ("non-homologous, excluding fragments", cross[~cross.frag], "tm"),
                       ("non-homologous, hand-curated family labels",
                        core[~core.same_family], "tm")]:
        r, p = mantel(g, x=x)
        res.append(dict(stratum=name, predictor=x, n=len(g), spearman=r, p_perm=p))
        print(f"{name}: n={len(g)} Spearman({x}, logOR)={r:.3f}, permutation p={p:.4f}")
    pd.DataFrame(res).to_csv(OUT / "structure_mantel.csv", index=False)

    # 非同源但結構像的對逐一列出，供人工檢查
    cross[cross.tm >= TM_HI].sort_values("z", ascending=False)[
        ["a", "b", "tm", "seqid80", "beta", "z", "fam_a", "fam_b"]
    ].to_csv(OUT / "nonhomologous_tm_hi.csv", index=False)


if __name__ == "__main__":
    main()
