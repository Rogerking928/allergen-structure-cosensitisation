"""聯合模型裡結構與序列兩個係數「差多少」，不是各自的 P 值落在 0.05 哪一邊。

外部審閱指出：TM P=0.023、序列 P=0.072 只是兩個 P 值分在 0.05 兩側，
不構成「結構比序列多帶資訊」的檢定（Gelman & Stern 2006）。這裡直接估 β_TM − β_序列，
信賴區間用節點自助法（同 discrim.py：抽成分、一對的權重＝兩端被抽中次數相乘），
差值的 P 不另算：逐對交換的排列忽略了對之間的相依，會太小。
"""
import json
import numpy as np
from scipy import stats
from common import ROOT
import robust
import discrim

OUT = ROOT / "out"
RNG = np.random.default_rng(20261008)
N = 5000


def design(g):
    X = np.column_stack([np.ones(len(g)), stats.rankdata(g.tm), stats.rankdata(g.seqid80)])   # 80-aa identity，與 robust.SEQ 一致
    X[:, 1:] = (X[:, 1:] - X[:, 1:].mean(0)) / X[:, 1:].std(0)
    return X


def fit(X, y, w=None):
    if w is not None:
        s = np.sqrt(w)
        X, y = X * s[:, None], y * s
    b, *_ = np.linalg.lstsq(X, y, rcond=None)
    return b[1], b[2]


def main():
    g = robust.homologous_core().reset_index(drop=True)
    X, y = design(g), g.beta.values
    b_tm, b_sq = fit(X, y)
    d_obs = b_tm - b_sq

    boots = []
    for _ in range(N):
        w = discrim.node_boot_weights(g, RNG).astype(float)
        if (w > 0).sum() < 30:
            continue
        t, s = fit(X, y, w)
        boots.append((t, s, t - s))
    boots = np.array(boots)
    lo, hi = np.percentile(boots, [2.5, 97.5], axis=0)


    res = dict(n=len(g), b_tm=b_tm, b_seq=b_sq, diff=d_obs,
               b_tm_ci=[lo[0], hi[0]], b_seq_ci=[lo[1], hi[1]], diff_ci=[lo[2], hi[2]],
               n_boot=len(boots))

    # 強共同致敏裡「兩個準則都不達標」那一格，拆成真正不同源與同源但兩者都不達標
    c_ = discrim.core_pairs()
    s_ = c_[c_.strong == 1]
    neither = s_[~s_.codex & ~s_.struct]
    res.update(strong=int(len(s_)), strong_nonhom=int((~s_.homologous).sum()),
               strong_neither=int(len(neither)),
               strong_neither_nonhom=int((~neither.homologous).sum()),
               strong_neither_hom=int(neither.homologous.sum()))
    (OUT / "coef_diff.json").write_text(json.dumps(res, indent=1))
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
