"""兩兩結構相似度（foldseek/TM-align）與序列相似度（mmseqs2）。

結構用 TM-align 模式全對全算，不走快速啟發式：成分只有一百多個，算得完，
而且本文的主張就建立在 TM-score 上，不該用近似值。
序列比對找不到的配對不是缺值，是「序列層級測不到相似性」，identity 記 0。
"""
import subprocess, itertools
import numpy as np
import pandas as pd
from common import ROOT

AFDB = ROOT / "data" / "afdb"
OUT = ROOT / "out"
TMP = ROOT / "data" / "tmp_search"


def run(cmd):
    print("$", " ".join(map(str, cmd)))
    subprocess.run(cmd, check=True)


def pairwise_identity(m):
    """每一對的局部比對 identity，以及 FAO/WHO 的 80-aa 窗判準。"""
    from Bio import Align
    from Bio.Align import substitution_matrices
    al = Align.PairwiseAligner(mode="local", open_gap_score=-11, extend_gap_score=-1)
    al.substitution_matrix = substitution_matrices.load("BLOSUM62")
    seqs = dict(m[["accession", "sequence"]].drop_duplicates("accession").itertuples(index=False))
    rows = []
    accs = sorted(seqs)
    for i, a in enumerate(accs):
        for b in accs[i + 1:]:
            aln = al.align(seqs[a], seqs[b])[0]
            qa, qb = aln[0], aln[1]
            match = [x == y and x != "-" for x, y in zip(qa, qb)]
            alnlen = sum(1 for x, y in zip(qa, qb) if x != "-" and y != "-")
            ident = 100 * sum(match) / alnlen if alnlen else 0.0
            # 80 殘基滑窗的最高 identity（FAO/WHO Codex 的 35% 判準用的就是這個）
            best80 = 0.0
            if len(match) >= 80:
                c = np.cumsum([0] + [int(x) for x in match])
                best80 = max(100 * (c[k + 80] - c[k]) / 80 for k in range(len(match) - 79))
            rows.append(dict(a=a, b=b, seqid=ident, alnlen=alnlen, seqid_best80=best80))
    return pd.DataFrame(rows)


def main():
    OUT.mkdir(exist_ok=True)
    TMP.mkdir(parents=True, exist_ok=True)
    m = pd.read_csv(ROOT / "out" / "component_uniprot.csv")
    m = m[m.status == "ok"].dropna(subset=["accession"])
    have = {p.stem for p in AFDB.glob("*.pdb")}
    m = m[m.accession.isin(have)]

    fs = OUT / "foldseek_hits.tsv"
    (TMP / "fs").mkdir(parents=True, exist_ok=True)
    # 16 執行緒同時配置 TMalign 緩衝區會在這台機器上 OOM（overcommit_memory=0），
    # 所以限 4 緒並設上限；157 個結構全對全本來就只要幾分鐘
    if not fs.exists():
        run(["foldseek", "easy-search", str(AFDB), str(AFDB), str(fs), str(TMP / "fs"),
             "--alignment-type", "1", "--exhaustive-search", "1", "-e", "inf",
             "--threads", "4", "--split-memory-limit", "2G",
             "--format-output", "query,target,fident,alntmscore,qtmscore,ttmscore,lddt,evalue",
             "--max-seqs", "4000"])
    f = pd.read_csv(fs, sep="\t", header=None,
                    names=["q", "t", "fident", "tm", "qtm", "ttm", "lddt", "evalue"])
    f["q"] = f.q.str.replace(".pdb", "", regex=False)
    f["t"] = f.t.str.replace(".pdb", "", regex=False)

    # 序列相似度用實際的兩兩比對算，不走 mmseqs：
    # 這裡只有一百多條序列，mmseqs 的預設配置在這台機器上反而會 OOM，
    # 而且交叉反應文獻慣用的判準（FAO/WHO：80 個胺基酸窗內 identity ≥ 35%）
    # 本來就需要真正的比對，不是偵測得到與否。
    s = pairwise_identity(m)
    s.to_csv(OUT / "pair_seqid.csv", index=False)

    # 每個無序配對留一列；TM-score 取兩個方向的平均（foldseek 的 tm 是對稱化後的值，
    # qtm/ttm 不對稱，另外保留最大值供敏感度分析）
    accs = sorted(m.accession.unique())
    rows = []
    fi = f.set_index(["q", "t"])
    si = s.set_index(["a", "b"])
    for a, b in itertools.combinations(accs, 2):
        tm = [fi.tm.get((x, y)) for x, y in [(a, b), (b, a)] if (x, y) in fi.index]
        tmax = [fi.qtm.get((x, y)) for x, y in [(a, b), (b, a)] if (x, y) in fi.index]
        sid = si.loc[(a, b)] if (a, b) in si.index else None
        rows.append(dict(a=a, b=b,
                         tm=pd.Series(tm).mean() if tm else 0.0,
                         tm_max=pd.Series(tmax).max() if tmax else 0.0,
                         # 以較長鏈正規化（兩個方向取小）。長度差很大時這才是保守的，
                         # 取大的話 21 殘基的片段會對上任何東西都 >0.9
                         tm_min=pd.Series(tmax).min() if tmax else 0.0,
                         seqid=float(sid.seqid) if sid is not None else 0.0,
                         seqid_best80=float(sid.seqid_best80) if sid is not None else 0.0,
                         alnlen=int(sid.alnlen) if sid is not None else 0))
    d = pd.DataFrame(rows)
    d.to_csv(OUT / "pair_similarity.csv", index=False)
    print(d.describe().to_string())
    print("結構像（TM≥0.5）但序列不像（80-aa 窗 <35%）:",
          int(((d.tm >= 0.5) & (d.seqid_best80 < 35)).sum()), "對")


if __name__ == "__main__":
    main()
