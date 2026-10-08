"""兩個過敏原的 AlphaFold 模型疊合（主鏈 Cα 軌跡），給 Figure 4 用。

為什麼不直接用主管線那次 foldseek 的結果：主管線跑的是 3Di 的快速比對（預設模式），
只拿 TM-score，不需要疊合矩陣。要畫圖就得有正確的旋轉平移，所以這裡用
`--alignment-type 1`（TM-align 模式）重跑這幾對，並把疊合後的 RMSD 印出來驗算——
3Di 模式算出來的對應關係疊起來 RMSD 11 Å，TM-align 模式 3.4 Å，不檢查會畫出錯的圖。

圖上標的 TM-score 一律用**主分析那個值**（以長鏈正規化），不是這裡重跑的值，
兩者定義不同，混用會讓正文與圖對不起來。
"""
import json
import subprocess
import sys
import numpy as np
import pandas as pd
from common import ROOT
import families

OUT = ROOT / "out"
AFDB = ROOT / "data" / "afdb"
TMP = ROOT / "data" / "tmp_overlay"

PAIRS = [("Cry_j_1", "Pla_a_2"), ("Ara_h_2", "Api_g_2")]


def ca(acc):
    p = AFDB / f"{acc}.pdb"
    return np.array([[float(l[30:38]), float(l[38:46]), float(l[46:54])]
                     for l in p.read_text().splitlines()
                     if l.startswith("ATOM") and l[12:16].strip() == "CA"])


def align(q_acc, t_acc):
    TMP.mkdir(parents=True, exist_ok=True)
    out = TMP / f"{q_acc}_{t_acc}.tsv"
    subprocess.run(["foldseek", "easy-search", str(AFDB / f"{q_acc}.pdb"),
                    str(AFDB / f"{t_acc}.pdb"), str(out), str(TMP / "work"),
                    "--format-output", "query,target,alntmscore,qstart,tstart,qaln,taln,u,t",
                    "-e", "1000", "--exhaustive-search", "1", "--alignment-type", "1",
                    "--threads", "2"],
                   check=True, capture_output=True)
    r = out.read_text().strip().split("\n")[0].split("\t")
    qs, ts, qa, ta = int(r[3]), int(r[4]), r[5], r[6]
    u = np.array([float(x) for x in r[7].split(",")]).reshape(3, 3)
    t = np.array([float(x) for x in r[8].split(",")])
    qi, ti, pq, pt = qs - 1, ts - 1, [], []
    for a, b in zip(qa, ta):
        if a != "-" and b != "-":
            pq.append(qi); pt.append(ti)
        qi += a != "-"
        ti += b != "-"
    return u, t, pq, pt


def main():
    m = families.table()
    m = m[(m.status == "ok") & m.accession.notna()].drop_duplicates("component")
    acc = m.set_index("component").accession
    p = pd.read_csv(OUT / "pairs_structure_cosens.csv")
    meta = {}
    for a, b in PAIRS:
        r = p[((p.a == a) & (p.b == b)) | ((p.a == b) & (p.b == a))].iloc[0]
        u, t, pq, pt = align(acc[a], acc[b])
        Q, T = ca(acc[a]), ca(acc[b])
        T2 = T @ u.T + t
        rmsd = float(np.sqrt(((Q[pq] - T2[pt]) ** 2).sum(1).mean()))
        if rmsd > 6:
            raise SystemExit(f"{a}-{b} 疊合 RMSD {rmsd:.1f} Å，太差，不要畫")
        np.savez(OUT / f"overlay_{a}_{b}.npz", q=Q, t=T2, u=u, tr=t, pq=pq, pt=pt,
                 acc_q=acc[a], acc_t=acc[b])
        meta[f"{a}|{b}"] = dict(tm=float(r.tm), seqid80=float(r.seqid80),
                                odds_ratio=float(np.exp(r.beta)), z=float(r.z),
                                aligned=len(pq), rmsd=rmsd,
                                homologous=bool(r.homologous))
        print(f"{a} - {b}: TM {r.tm:.2f}, identity {r.seqid80:.0f}%, "
              f"OR {np.exp(r.beta):.1f}, aligned {len(pq)}, RMSD {rmsd:.2f} Å")
    (OUT / "overlay.json").write_text(json.dumps(meta, indent=1))


if __name__ == "__main__":
    main()
