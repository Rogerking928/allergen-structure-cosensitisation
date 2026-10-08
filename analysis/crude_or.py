"""未校正的兩兩共同致敏 log OR（2×2，Haldane–Anscombe 修正）。

只為了對照：用來顯示不校正個人致敏廣度時，整張表會偏成一片正相關。
矩陣相乘一次算完，不逐對配適。
"""
import itertools
import numpy as np
import pandas as pd
from common import load, ROOT
from cosens import POS, MIN_POS, MIN_CELL, strat_table


def main():
    d, meta, comps = load()
    out = []
    for chip in d.Chip_Type.unique():
        y, keep = strat_table(d, comps, chip)
        npos = y.sum(axis=0)
        use = [c for c in keep if npos[c] >= MIN_POS]
        a = np.nan_to_num(y[use].values)
        obs = (~np.isnan(y[use].values)).astype(float)
        n11 = a.T @ a
        n10 = a.T @ (obs - a)
        n01 = (obs - a).T @ a
        n00 = (obs - a).T @ (obs - a)
        for i, j in itertools.combinations(range(len(use)), 2):
            if n11[i, j] < MIN_CELL:
                continue
            or_ = ((n11[i, j] + .5) * (n00[i, j] + .5)) / ((n10[i, j] + .5) * (n01[i, j] + .5))
            se = np.sqrt(sum(1 / (x[i, j] + .5) for x in (n11, n10, n01, n00)))
            out.append(dict(chip=chip, a=use[i], b=use[j], beta_crude=np.log(or_), se_crude=se,
                            n11=int(n11[i, j])))
    r = pd.DataFrame(out)
    r.to_csv(ROOT / "out" / "cosens_crude.csv", index=False)
    g = r.assign(w=1 / r.se_crude ** 2).groupby(["a", "b"])
    pooled = pd.DataFrame({
        "beta_crude": g.apply(lambda x: np.average(x.beta_crude, weights=x.w), include_groups=False),
        "se_crude": 1 / np.sqrt(g.w.sum()),
    }).reset_index()
    pooled.to_csv(ROOT / "out" / "cosens_crude_pooled.csv", index=False)
    print(pooled.beta_crude.describe().to_string())


if __name__ == "__main__":
    main()
