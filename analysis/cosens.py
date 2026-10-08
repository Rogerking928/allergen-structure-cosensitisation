"""逐對成分的共同致敏強度，校正個人整體致敏廣度。

為什麼要校正廣度：很會過敏的人對什麼都陽性，不校正的話任兩個成分都「常一起陽性」，
整張表會變成一片正相關，跟結構毫無關係也測得到。所以每一對都用
    logit P(Y_i=1) = a + b*Y_j + c*breadth
的 b 當結果，breadth 是該人在同一張晶片上、扣掉 i 與 j 之後的陽性成分數。

晶片型別分層（ISAC v1 / ISAC v2 / ALEX 的動態範圍與成分套組不同），
事後用逆變異數加權合併；只在單一型別測到的成分對照樣保留並標記。
"""
import itertools
import numpy as np
import pandas as pd
import statsmodels.api as sm
from common import load, ROOT

POS = 0.3          # ISU-E / kUA/L 的慣用陽性閾值
MIN_POS = 15       # 每個成分在該層至少要有這麼多陽性才估計
MIN_CELL = 5       # 共陽性格至少要有這麼多人


def strat_table(d, comps, chip):
    g = d[d.Chip_Type == chip]
    x = g[comps].replace(-1, np.nan)
    keep = [c for c in comps if x[c].notna().sum() >= 0.8 * len(g)]
    y = (x[keep] >= POS).astype(float).where(x[keep].notna())
    return y.dropna(axis=0, how="all"), keep


def fit_pair(yi, yj, breadth):
    X = sm.add_constant(np.column_stack([yj, breadth]))
    try:
        r = sm.Logit(yi, X).fit(disp=0)
        b, se = r.params[1], r.bse[1]
        if not np.isfinite(b) or not np.isfinite(se) or se > 5:
            raise ValueError
        return b, se, "mle"
    except Exception:
        r = sm.Logit(yi, X).fit_regularized(disp=0, alpha=1.0, L1_wt=0.0)
        return r.params[1], np.nan, "ridge"


def main():
    d, meta, comps = load()
    out = []
    for chip in d.Chip_Type.unique():
        y, keep = strat_table(d, comps, chip)
        tot = y.sum(axis=1).values
        npos = y.sum(axis=0)
        use = [c for c in keep if npos[c] >= MIN_POS]
        print(f"{chip}: n={len(y)} 成分 {len(keep)} → 可估 {len(use)}")
        arr = y[use].values
        idx = {c: k for k, c in enumerate(use)}
        for a, b in itertools.combinations(use, 2):
            ia, ib = idx[a], idx[b]
            ya, yb = arr[:, ia], arr[:, ib]
            m = ~(np.isnan(ya) | np.isnan(yb))
            if m.sum() < 50 or ((ya[m] == 1) & (yb[m] == 1)).sum() < MIN_CELL:
                continue
            breadth = tot[m] - ya[m] - yb[m]
            beta, se, how = fit_pair(ya[m], yb[m], breadth)
            out.append(dict(chip=chip, a=a, b=b, n=int(m.sum()),
                            n_a=int(ya[m].sum()), n_b=int(yb[m].sum()),
                            n_ab=int(((ya[m] == 1) & (yb[m] == 1)).sum()),
                            beta=beta, se=se, method=how))
    r = pd.DataFrame(out)
    r.to_csv(ROOT / "out" / "cosens_by_chip.csv", index=False)

    # 逆變異數合併（只用有標準誤的 MLE 估計）
    ok = r[(r.method == "mle") & r.se.notna()].copy()
    ok["w"] = 1 / ok.se ** 2
    g = ok.groupby(["a", "b"])
    pooled = pd.DataFrame({
        "beta": g.apply(lambda x: np.average(x.beta, weights=x.w), include_groups=False),
        "se": 1 / np.sqrt(g.w.sum()),
        "k_chip": g.size(),
        "n_total": g.n.sum(),
    }).reset_index()
    pooled["z"] = pooled.beta / pooled.se
    pooled.to_csv(ROOT / "out" / "cosens_pooled.csv", index=False)
    print(pooled.describe().to_string())


if __name__ == "__main__":
    main()
