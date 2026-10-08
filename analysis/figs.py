"""圖表本的圖（700 dpi PNG，無 PDF）。數字一律讀 out/*.csv／*.json，不在這裡重算。

2026-10-07 改版（Article 規格，四張主圖）：
  Figure 1  資料怎麼來、為什麼要校正致敏廣度（原本散在補充）
  Figure 2  相似度與共同致敏的關係（原 Figure 1，加上分帶的劑量反應）
  Figure 3  兩個篩選門檻：序列準則與結構門檻各抓到什麼
  Figure 4  生物學對照：網路、家族×家族、留一同源群
補充：S1 敏感度森林圖、S2 AlphaFold 模型信心、S3 排列檢定虛無分布、S4 三個平台。
"""
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, Rectangle
from matplotlib.lines import Line2D
from matplotlib.colors import TwoSlopeNorm, LogNorm
import figstyle as F
from common import load, ROOT
import families

OUT = ROOT / "out"

# 家族顏色：四張圖共用，讀者在 Figure 4 學到的對應不能改
FAM_COLORS = {
    "PR-10 (Bet v 1-like)": F.NAVY,
    "nsLTP": F.RED,
    "2S albumin (prolamin)": F.AMBER,
    "Profilin": F.GREEN,
    "Tropomyosin": F.PURPLE,
    "Lipocalin": F.TEAL,
    "11S globulin (legumin)": F.BROWN,
    "7S globulin (vicilin)": F.PLUM,
}
OTHER = "#b7c0c8"


def summ():
    return json.load(open(OUT / "summary.json"))


def load_pairs():
    adj = pd.read_csv(OUT / "cosens_pooled.csv")
    cru = pd.read_csv(OUT / "cosens_crude_pooled.csv")
    d = adj.merge(cru, on=["a", "b"], how="inner")
    fam = families.table()
    fam = fam[fam.status == "ok"].drop_duplicates("component").set_index("component")
    for s in "ab":
        d[f"fam_{s}"] = d[s].map(fam.family)
        d[f"org_{s}"] = d[s].map(fam.organism)
    d["same_family"] = (d.fam_a == d.fam_b) & d.fam_a.notna() & (d.fam_a != "Other / unclassified")
    d["same_species"] = (d.org_a == d.org_b) & d.org_a.notna()
    return d


def core():
    import discrim
    return discrim.core_pairs()


def _pairs():
    import structure_cosens as S
    p = pd.read_csv(OUT / "pairs_structure_cosens.csv")
    c = p[(~p.same_species) & (~p.subchain) & (p.len_a >= S.MIN_LEN) & (p.len_b >= S.MIN_LEN)
          & p.annotated]
    ss = p[p.same_species & ~p.subchain]
    return p, c, ss


def overlap_check(fig, skip=()):
    """文字框兩兩不可重疊，用 bbox 檢查而不是用眼睛看（skip：示意圖那種刻意貼齊的）。"""
    fig.canvas.draw()
    rr = fig.canvas.get_renderer()
    boxes = []
    for ax in fig.axes:
        if ax in skip:
            continue
        x0, x1 = sorted(ax.get_xlim()); y0, y1 = sorted(ax.get_ylim())
        items = [t for t in ax.get_xticklabels() if x0 <= t.get_position()[0] <= x1]
        items += [t for t in ax.get_yticklabels() if y0 <= t.get_position()[1] <= y1]
        items += list(ax.texts) + [ax.xaxis.label, ax.yaxis.label]
        leg = ax.get_legend()
        if leg:
            items.append(leg)
        for t in items:
            if hasattr(t, "get_text") and not t.get_text():
                continue
            if not t.get_visible():
                continue
            bb = t.get_window_extent(rr)
            if (bb.x0, bb.y0, bb.x1, bb.y1) == (0.0, 0.0, 1.0, 1.0):
                raise SystemExit(f"文字沒有被排版（多半是畫到了軸範圍外）：{t.get_text()!r}")
            boxes.append((t, bb))
    bad = []
    for i in range(len(boxes)):
        for j in range(i + 1, len(boxes)):
            a, b = boxes[i][1], boxes[j][1]
            if a.overlaps(b) and a.width > 0 and b.width > 0:
                bad.append((getattr(boxes[i][0], "get_text", lambda: "legend")(),
                            getattr(boxes[j][0], "get_text", lambda: "legend")()))
    if bad:
        raise SystemExit(f"文字重疊：{bad[:5]}")


def data_overlap_check(fig, skip=(), pad=2.0):
    """文字（註記、圖例的字）不可壓到資料：線（含 axvline／axhline 與誤差線）、散點、長條、hexbin。

    overlap_check 只管字與字；這裡管字與資料。用畫出來的像素座標算，不靠眼睛。"""
    from matplotlib.collections import PathCollection, PolyCollection, LineCollection
    from matplotlib.patches import Rectangle as _R
    fig.canvas.draw()
    rr = fig.canvas.get_renderer()
    bad = []

    def segs_hit(xy, bb):
        for p0, p1 in zip(xy[:-1], xy[1:]):
            t = np.linspace(0, 1, 60)[:, None]
            pts = p0 + t * (p1 - p0)
            if ((pts[:, 0] > bb.x0) & (pts[:, 0] < bb.x1) &
                    (pts[:, 1] > bb.y0) & (pts[:, 1] < bb.y1)).any():
                return True
        return False

    for ax in fig.axes:
        if ax in skip:
            continue
        texts = [t for t in ax.texts if t.get_text().strip() and t.get_visible()
                 and t.get_gid() != "inside-bar"]                 # 刻意寫在長條內的字
        leg = ax.get_legend()
        if leg is not None:
            texts += [t for t in leg.get_texts() if t.get_text().strip()]
            if leg.get_title().get_text().strip():
                texts.append(leg.get_title())
        for t in texts:
            bb = t.get_window_extent(rr).padded(pad)
            name = t.get_text().split("\n")[0][:40]
            for ln in ax.lines:
                if not ln.get_visible():
                    continue
                xy = ln.get_transform().transform(ln.get_xydata())
                xy = xy[np.isfinite(xy).all(1)]
                if ln.get_linestyle() in ("None", "", " "):      # 只有標記的線
                    hit = ((xy[:, 0] > bb.x0) & (xy[:, 0] < bb.x1) &
                           (xy[:, 1] > bb.y0) & (xy[:, 1] < bb.y1)).any()
                else:
                    hit = len(xy) > 1 and segs_hit(xy, bb)
                if hit:
                    bad.append((name, "line"))
            for col in ax.collections:
                if not col.get_visible():
                    continue
                if isinstance(col, LineCollection):
                    for seg in col.get_segments():
                        if len(seg) > 1 and segs_hit(col.get_transform().transform(seg), bb):
                            bad.append((name, "errorbar")); break
                elif isinstance(col, (PathCollection, PolyCollection)):
                    off = col.get_offsets()
                    if len(off) == 0:
                        continue
                    if isinstance(col, PolyCollection):           # hexbin：位移在資料座標
                        xy = ax.transData.transform(off)
                        v = col.get_paths()[0].vertices
                        d0 = ax.transData.transform(off[:1])[0]
                        r = 0.5 * np.ptp(ax.transData.transform(off[0] + v) - d0, axis=0).max()
                    else:
                        xy = col.get_offset_transform().transform(off)
                        sz = col.get_sizes()
                        r = np.sqrt(sz.max() if len(sz) else 1) / 2 * fig.dpi / 72
                    near = ((xy[:, 0] > bb.x0 - r) & (xy[:, 0] < bb.x1 + r) &
                            (xy[:, 1] > bb.y0 - r) & (xy[:, 1] < bb.y1 + r))
                    if near.any():
                        bad.append((name, f"points ({near.sum()})"))
            for pa in ax.patches:
                if isinstance(pa, _R) and pa.get_visible() and pa.get_width() != 0 \
                        and pa.get_window_extent(rr).overlaps(bb):
                    bad.append((name, "bar"))
    # 相鄰兩格的軸標題也不可互碰（S1 的 A、B 橫軸標題曾經黏在一起）
    labs = [t for a_ in fig.axes for t in (a_.xaxis.label, a_.yaxis.label) if t.get_text()]
    for i_ in range(len(labs)):
        for j_ in range(i_ + 1, len(labs)):
            if labs[i_].get_window_extent(rr).padded(pad).overlaps(
                    labs[j_].get_window_extent(rr).padded(pad)):
                bad.append((labs[i_].get_text()[:30], labs[j_].get_text()[:30]))
    if bad:
        raise SystemExit(f"文字壓到資料：{sorted(set(bad))[:8]}")


def or_axis(ax, lo=-6.2, hi=9.5):
    """y 軸用 log OR 作圖、但刻度標勝算比：臨床讀者腦中是 OR，不是 log OR。"""
    ticks = [0.01, 0.1, 1, 10, 100, 1000]
    ax.set_ylim(lo, hi)
    ax.set_yticks([np.log(t) for t in ticks if lo <= np.log(t) <= hi])
    ax.set_yticklabels([("%g" % t) for t in ticks if lo <= np.log(t) <= hi])


# ============================================================ Figure 1（方法）
def _box(ax, x, y, w, h, text, fc="white", ec=F.NAVY, fs=7.2):
    r = ax.add_patch(Rectangle((x, y), w, h, fc=fc, ec=ec, lw=0.9, zorder=2))
    t = ax.annotate(text, xy=(x + w / 2, y + h / 2), ha="center", va="center", fontsize=fs,
                    zorder=3)
    return r, t


def box_text_inside(fig, boxes, pad=3.0):
    """示意圖的字必須整個在框內、離框線至少 pad 像素（overlap_check 跳過示意圖，這裡補上）。"""
    fig.canvas.draw()
    rr = fig.canvas.get_renderer()
    for r, t in boxes:
        bb, tb = r.get_window_extent(rr), t.get_window_extent(rr)
        if (tb.x0 < bb.x0 + pad or tb.x1 > bb.x1 - pad or
                tb.y0 < bb.y0 + pad or tb.y1 > bb.y1 - pad):
            raise SystemExit(f"示意圖的字壓到框線：{t.get_text()!r} 框 {bb} 字 {tb}")


def _arrow(ax, p0, p1):
    ax.add_patch(FancyArrowPatch(p0, p1, arrowstyle="-|>", mutation_scale=7,
                                 lw=0.9, color=F.GREY, shrinkA=0, shrinkB=0, zorder=1))


def figure1():
    s = summ()
    d, meta, comps = load()
    fig = plt.figure(figsize=(F.WIDTH_IN, 5.0))
    gs = fig.add_gridspec(2, 2, height_ratios=[0.92, 1.0], hspace=0.38, wspace=0.3)

    # --- A 示意圖：箭頭一律停在框外（量過 bbox，不靠目測）
    ax = fig.add_subplot(gs[0, 0])
    ax.set_xlim(0, 10); ax.set_ylim(1.55, 9.58); ax.axis("off")
    steps = [(f"{s['n_patients']:,} patients, 3 platforms\n(ISAC v1, ISAC v2, ALEX2)", 8.0),
             (f"{s['n_component_columns']} IgE columns → {s['n_mapped']} single-protein\n"
              f"analytes (WHO/IUIS → UniProt)", 5.9),
             ("AlphaFold models → TM-align\nand 80-aa window identity", 3.8),
             ("Pair-level odds ratio,\nadjusted for sensitisation breadth", 1.7)]
    boxes = [_box(ax, 0.4, y, 9.2, 1.5, text) for text, y in steps]
    for _, y in steps[:-1]:
        _arrow(ax, (5.0, y - 0.04), (5.0, y - 0.55))
    F.panel(ax, "A", dx=-6, dy=-4)

    # --- B 漏斗（S2 表的數字）：線性、左對齊、有 x 軸——長條的長度就是數量，
    #     最大的一刀（缺 AlphaFold 模型）要看得出來；對數長度沒有軸會被讀成線性
    ax = fig.add_subplot(gs[0, 1])
    rows = [("pairs estimated", s["pairs_pooled"]),
            ("with a model for both", s["pairs_with_structure"]),
            ("analysed pairs", s["pairs_core"]),
            ("homologous", s["pairs_homologous"])]
    excl = [f"\u2212 {s['pairs_pooled'] - s['pairs_with_structure']:,} no model for one member",
            f"\u2212 {s['pairs_with_structure'] - s['pairs_core']:,} same species, sub-chain, "
            "< 50 aa, unannotated",
            f"\u2212 {s['pairs_nonhomologous']:,} non-homologous"]
    xmax = 20000
    for i, (lab, n) in enumerate(rows):
        y = 3 - i
        ax.barh(y, n, height=0.42, color=F.NAVY, alpha=0.85, lw=0)
        if n > 0.6 * xmax:     # 第一根太長，字寫在條內
            ax.annotate(f"{n:,} {lab}", xy=(300, y), ha="left", va="center", fontsize=6.8,
                        color="white", gid="inside-bar")
        else:
            ax.annotate(f"{n:,} {lab}", xy=(n + 250, y), ha="left", va="center",
                        fontsize=6.8)
        if i < 3:
            ax.annotate(excl[i], xy=(250, y - 0.5), ha="left", va="center",
                        fontsize=6.2, color=F.GREY)
    ax.set_yticks([])
    ax.set_xlim(0, xmax); ax.set_ylim(-0.4, 3.4)
    ax.set_xticks([0, 5000, 10000, 15000, 20000])
    ax.set_xticklabels(["0", "5,000", "10,000", "15,000", "20,000"])
    ax.set_xlabel("Allergen pairs")
    F.bare(ax)
    ax.spines["left"].set_visible(False)
    F.panel(ax, "B", dx=-14, dy=-4)

    # --- C 致敏廣度
    ax = fig.add_subplot(gs[1, 0])
    for chip in F.CHIP_ORDER:
        g = d[d.Chip_Type == chip]
        x = g[comps].replace(-1, np.nan)
        keep = [c for c in comps if x[c].notna().sum() >= 0.8 * len(g)]
        breadth = 100 * (x[keep] >= 0.3).sum(axis=1) / len(keep)
        bins = np.arange(0, 61, 2)
        # 超過 60% 的病人併進最後一格（每平台各 1 人，原本被丟掉、合計少 1 人）
        h, e = np.histogram(breadth.clip(upper=59.9), bins=bins)
        assert h.sum() == len(breadth)
        ax.step(e[:-1], 100 * h / len(breadth), where="post", color=F.CHIP[chip], lw=1.3,
                label=f"{F.CHIP_LABEL[chip]} ({len(g):,})")
    F.bare(ax)
    ax.set_xlabel("Components positive per patient (%)")
    ax.set_ylabel("Patients (%)")
    ax.legend(frameon=False, loc="upper right", fontsize=7.5)
    F.panel(ax, "C", dx=-34)

    # --- D 校正前後
    ax = fig.add_subplot(gs[1, 1])
    p = load_pairs()
    ax.scatter(p.beta_crude, p.beta, s=3, lw=0, alpha=0.25, color=F.METHOD["main"],
               rasterized=True)
    lim = [-4, 8.5]
    ax.plot(lim, lim, lw=0.8, color=F.METHOD["other"], ls=(0, (4, 3)))
    ax.axhline(0, lw=0.8, color=F.METHOD["other"])
    ax.set_xlim(lim); ax.set_ylim(lim)
    F.bare(ax)
    ax.set_xlabel("Unadjusted log OR")
    ax.set_ylabel("Breadth-adjusted log OR")
    ax.annotate(f"median {p.beta_crude.median():.2f} → {p.beta.median():.2f}\n"
                f"{len(p):,} pairs", xy=(0.04, 0.90), xycoords="axes fraction", fontsize=7.5)
    F.panel(ax, "D", dx=-30)
    overlap_check(fig, skip=(fig.axes[0],))
    data_overlap_check(fig, skip=(fig.axes[0],))
    box_text_inside(fig, boxes)
    F.save(fig, "Figure_1")


# ============================================================ Figure 2（主結果）
def figure2():
    s = summ()
    mt = {r["stratum"]: r for r in s["mantel"]}
    b = json.load(open(OUT / "bands.json"))
    bands = pd.read_csv(OUT / "tm_bands_ci.csv")
    p, c, ss = _pairs()
    hom, non = c[c.homologous], c[~c.homologous]

    fig = plt.figure(figsize=(F.WIDTH_IN, 4.8))
    gs = fig.add_gridspec(2, 6, hspace=0.55, wspace=0.78, right=0.80)

    ax = fig.add_subplot(gs[0, 0:2])
    rng = np.random.default_rng(0)
    groups = [("Same\nspecies", ss), ("Homo-\nlogous", hom), ("Non-\nhomologous", non)]
    for i, (name, g) in enumerate(groups):
        y = g.beta.values
        ax.scatter(rng.normal(i, 0.08, len(y)), y, s=2, lw=0,
                   alpha=0.18 if len(y) > 2000 else 0.35, color=F.METHOD["other"], rasterized=True)
        q = np.percentile(y, [25, 50, 75])
        ax.plot([i - 0.3, i + 0.3], [q[1], q[1]], lw=1.6, color=F.INK, zorder=3)
        ax.plot([i, i], [q[0], q[2]], lw=1.1, color=F.INK, zorder=3)
        ax.annotate(f"{len(y):,}", xy=(i, -5.9), fontsize=6.6, ha="center", color=F.GREY)
    ax.axhline(0, lw=0.7, color=F.METHOD["other"], ls=(0, (4, 3)))
    ax.set_xticks(range(3)); ax.set_xticklabels([g[0] for g in groups], fontsize=6.0)
    ax.set_xlim(-0.55, 2.55)
    or_axis(ax)
    ax.set_ylabel("Co-sensitisation\n(odds ratio)")
    F.bare(ax); F.panel(ax, "A", dx=-38)

    # B 同源對 vs TM，點大小反映精確度（1/SE）
    ax = fig.add_subplot(gs[0, 2:])
    hi = hom.seqid80 >= 35
    sz = np.clip(6 / hom.se, 2, 26)
    ax.scatter(hom.tm[~hi], hom.beta[~hi], s=sz[~hi], lw=0.5, facecolor="white",
               edgecolor=F.NAVY, label="80-aa identity < 35%")
    ax.scatter(hom.tm[hi], hom.beta[hi], s=sz[hi], lw=0, color=F.NAVY, alpha=0.75,
               label="80-aa identity ≥ 35%")
    for r in s["tm_bins"]:
        lo, up = [float(x) for x in r["tm_bin"].strip("[)").split(",")]
        ax.plot([lo, min(up, 1.0)], [r["median_logOR"]] * 2, color=F.RED, lw=1.8, zorder=4)
    ax.axhline(0, lw=0.7, color=F.METHOD["other"], ls=(0, (4, 3)))
    ax.set_xlim(0.15, 1.0)
    or_axis(ax)
    ax.set_xlabel("Structural similarity (TM-score)")
    r = mt["homologous"]
    # 圖例放右下（同源對在高 TM 處的 OR 都在 1 以上，那裡是空的）；左上只留一行說明
    ax.legend(frameon=False, loc="lower right", fontsize=7, handletextpad=0.2, borderaxespad=0.1)
    ax.annotate(f"Homologous pairs, n = {r['n']};  ρ = {r['spearman']:.2f}", xy=(0.02, 0.97),
                xycoords="axes fraction", ha="left", va="top", fontsize=7.2)
    ax.annotate("point area ∝ 1/SE", xy=(0.02, 0.88), xycoords="axes fraction",
                ha="left", va="top", fontsize=6.6, color=F.GREY)
    F.bare(ax); F.panel(ax, "B", dx=-30)

    # C 非同源（x 軸只畫到 0.55：沒有任何一對超過 0.47，留白會讓人以為資料被截掉）
    ax = fig.add_subplot(gs[1, 0:3])
    ax.hexbin(non.tm, non.beta, gridsize=(26, 20), extent=(0.15, 0.55, -6.2, 9.5),
              cmap="Greys", mincnt=1, lw=0.1,
              norm=LogNorm(vmin=0.35))   # 對數色階的下限壓到 1 以下：count=1 原本是白色、等於沒畫
    ax.axhline(0, lw=0.7, color=F.METHOD["other"], ls=(0, (4, 3)))
    ax.axvline(0.5, lw=0.8, color=F.RED, ls=(0, (2, 2)))
    ax.set_xlim(0.15, 0.56)
    or_axis(ax)
    ax.set_xlabel("Structural similarity (TM-score)")
    ax.set_ylabel("Odds ratio")
    r = mt["non-homologous"]
    ax.annotate(f"Non-homologous pairs, n = {r['n']:,}\nρ = {r['spearman']:.3f} "
                f"(P = {r['p_perm']:.2f})\nhighest TM-score reached "
                f"{s['nonhomologous_tm_max']:.2f}",
                xy=(0.03, 0.96), xycoords="axes fraction", ha="left", va="top", fontsize=7.0)
    F.bare(ax); F.panel(ax, "C", dx=-30)

    # D 分帶劑量反應
    ax = fig.add_subplot(gs[1, 3:])
    y = np.arange(len(bands))[::-1]
    ax.errorbar(bands.median_or, y, xerr=[bands.median_or - bands.lo, bands.hi - bands.median_or],
                fmt="o", ms=3.4, lw=1.1, color=F.NAVY, capsize=0)
    ax.set_xscale("log")
    ax.set_xticks([1, 10, 100]); ax.set_xticklabels(["1", "10", "100"])
    ax.axvline(1, lw=0.7, color=F.METHOD["other"], ls=(0, (4, 3)))
    ax.yaxis.tick_right()
    ax.set_yticks(y); ax.set_yticklabels([f"{b} ({n})" for b, n in zip(bands.band, bands.n)],
                                         fontsize=6.8)
    ax.set_ylim(-0.7, len(bands) - 0.3)
    ax.set_xlabel("Median odds ratio (95% CI)", fontsize=8)
    ax.annotate(f"trend P {b['trend_p_homologous']:.4f}".replace("0.0001", "< 0.001"),
                xy=(0.97, 0.97), xycoords="axes fraction", ha="right", va="top",
                fontsize=6.8, color=F.GREY)
    F.bare(ax, keep=("bottom", "right")); F.panel(ax, "D", dx=-10)
    overlap_check(fig)
    data_overlap_check(fig)
    F.save(fig, "Figure_2")


# ============================================================ Figure 3（兩個門檻）
def figure3():
    """兩個門檻各抓到什麼。

    **敘事順序**：A 四象限 → B 每格的量 → C 具名配對 → D 判別力。
    判別力那格寫「兩者的 AUC 相當」與聯合模型係數差（分不出來）。
    序列用 80-aa 後，同源對內序列 ρ 0.82 ≥ 結構 0.79；結構的價值在 116 對與分半驗證，不在係數。
    """
    import structure_cosens as S
    dd = json.load(open(OUT / "discrimination.json"))
    sm = summ()
    rb = {(r["analysis"], r["term"]): r for r in sm["robust"]}
    c = core()
    q = pd.read_csv(OUT / "quadrant_struct_only.csv")
    top = q[q.z > 3].head(14)

    fig = plt.figure(figsize=(F.WIDTH_IN, 5.6))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.05, 1.0], hspace=0.62, wspace=0.52)

    # --- A 二維平面
    ax = fig.add_subplot(gs[0, 0])
    o = c.sort_values("beta")
    norm = TwoSlopeNorm(vmin=-3, vcenter=0, vmax=6)
    sc = ax.scatter(o.seqid80, o.tm, c=o.beta, cmap="RdYlBu_r", norm=norm, s=6, lw=0,
                    rasterized=True)
    named = set(zip(top.a, top.b))
    mk = c[[(x, y) in named or (y, x) in named for x, y in zip(c.a, c.b)]]
    ax.scatter(mk.seqid80, mk.tm, s=26, facecolor="none", edgecolor=F.INK, lw=0.7,
               label="shown in C")
    # 縱線只畫在資料範圍內：上下留白是給四格標籤用的，線穿過去就會壓到字
    ax.plot([S.SEQ_HI, S.SEQ_HI], [0.08, 1.04], lw=0.9, color=F.INK, ls=(0, (3, 2)))
    ax.axhline(S.TM_HI, lw=0.9, color=F.INK, ls=(0, (3, 2)))
    # 四格的標註一律放角落，墊白底：這個版面沒有真正的空白處，硬找位置會壓到點
    cells = {(True, False): (0.015, 0.985), (True, True): (0.985, 0.985),
             (False, False): (0.015, 0.015), (False, True): (0.985, 0.015)}
    for r in dd["reclassification"]:
        x, y = cells[(r["structural"], r["codex"])]
        ax.annotate(f"n = {r['n']:,}\nOR {r['median_or']:.1f}", xy=(x, y),
                    xycoords="axes fraction", ha="right" if x > 0.5 else "left",
                    va="top" if y > 0.5 else "bottom", fontsize=6.4, zorder=6,
                    bbox=dict(fc="white", alpha=0.78, lw=0, pad=1.4))
    ax.set_xlabel("Best 80-aa window identity (%)")
    ax.set_ylabel("Structural similarity (TM-score)")
    ax.set_xlim(-7, 103); ax.set_ylim(-0.2, 1.33)
    ax.set_yticks([0.2, 0.4, 0.6, 0.8, 1.0])
    ax.set_xticks([0, 20, 40, 60, 80, 100])
    ax.set_xticklabels(["none", "20", "40", "60", "80", "100"])
    ax.legend(frameon=False, loc="lower center", fontsize=6.2, handletextpad=0.2,
              borderaxespad=0.15, bbox_to_anchor=(0.5, -0.42))
    cb = fig.colorbar(sc, ax=ax, pad=0.02, fraction=0.045)
    cb.set_label("Adjusted log odds ratio", fontsize=7)
    cb.ax.tick_params(labelsize=6.8)
    cb.outline.set_linewidth(0.5)
    F.bare(ax); F.panel(ax, "A", dx=-36, dy=10)

    # --- B 每一格的量（本圖的重心）
    ax = fig.add_subplot(gs[0, 1])
    labs = {(True, True): "Both criteria", (True, False): "Structural\nthreshold only",
            (False, True): "Sequence\ncriterion only", (False, False): "Neither"}
    order = [(True, True), (True, False), (False, True), (False, False)]
    rows = sorted(dd["reclassification"],
                  key=lambda r: order.index((r["structural"], r["codex"])))
    y = np.arange(len(rows))[::-1]
    for yy, r in zip(y, rows):
        ax.barh(yy, r["n"], color=OTHER, height=0.62, lw=0)
        ax.barh(yy, r["strong"], color=F.RED, height=0.62, lw=0)
        ax.annotate(f"{r['strong']:,}/{r['n']:,} ({r['pct_strong']:.0f}%)",
                    xy=(r["n"] * 1.3, yy), va="center", fontsize=6.8)
    ax.set_xscale("log"); ax.set_xlim(1, 6e4)
    ax.set_yticks(y)
    ax.set_yticklabels([labs[(r["structural"], r["codex"])] for r in rows], fontsize=7.0)
    ax.set_ylim(-0.75, len(rows) - 0.25)
    ax.set_xlabel("Pairs (log scale)")
    ax.legend(handles=[Line2D([], [], marker="s", ls="", ms=6, color=F.RED,
                              label="strongly co-sensitised (z > 3)"),
                       Line2D([], [], marker="s", ls="", ms=6, color=OTHER, label="other")],
              frameon=False, loc="upper center", fontsize=6.6, borderaxespad=0.2,
              bbox_to_anchor=(0.45, -0.30), ncol=2, handletextpad=0.3)
    F.bare(ax); F.panel(ax, "B", dx=-60)

    # --- C 具名配對，右側附 TM 與序列相似度
    ax = fig.add_subplot(gs[1, 0])
    t = top.iloc[::-1]
    y = np.arange(len(t))
    orr = np.exp(t.beta)
    lo, hi = np.exp(t.beta - 1.96 * t.se), np.exp(t.beta + 1.96 * t.se)
    ax.errorbar(orr, y, xerr=[orr - lo, hi - orr], fmt="o", ms=3.0, lw=0.9,
                color=F.NAVY, capsize=0)
    ax.set_xscale("log"); ax.set_xlim(1.2, 900)
    ax.axvline(1, lw=0.7, color=F.METHOD["other"], ls=(0, (4, 3)))
    ax.set_yticks(y)
    ax.set_yticklabels([f"{a.replace('_RUO', '').replace('_', ' ')} \u2013 "
                        f"{b.replace('_RUO', '').replace('_', ' ')}"
                        for a, b in zip(t.a, t.b)], fontsize=6.4)
    ax.set_ylim(-1.3, len(t) + 0.5)
    ax.set_xlabel("Odds ratio (95% CI)")
    for yy, r in zip(y, t.itertuples()):
        ax.annotate(f"{r.tm:.2f}", xy=(1.04, yy), xycoords=("axes fraction", "data"),
                    ha="right", va="center", fontsize=6.0, color=F.GREY)
        ax.annotate(f"{r.seqid80:.0f}", xy=(1.17, yy), xycoords=("axes fraction", "data"),
                    ha="right", va="center", fontsize=6.0, color=F.GREY)
    ax.annotate("TM", xy=(1.04, len(t) - 0.25), xycoords=("axes fraction", "data"),
                ha="right", va="bottom", fontsize=6.0, color=F.GREY)
    ax.annotate("id %", xy=(1.17, len(t) - 0.25), xycoords=("axes fraction", "data"),
                ha="right", va="bottom", fontsize=6.0, color=F.GREY)
    F.bare(ax); F.panel(ax, "C", dx=-82)

    # --- D 判別力：兩條 ROC 相當，增量在聯合模型
    ax = fig.add_subplot(gs[1, 1])
    auc = {r["predictor"]: r for r in dd["auc_all"]}
    for col, colr, lab in [("tm", F.NAVY, "TM-score"), ("seqid80", F.AMBER, "80-aa identity")]:
        r = pd.read_csv(OUT / f"roc_{col}.csv")
        a = auc[col]
        ax.plot(r.fpr, r.tpr, lw=1.3, color=colr,
                label=f"{lab} {a['auc']:.2f} ({a['lo']:.2f}\u2013{a['hi']:.2f})")
    pts_h = []
    for pt, mk_, lab in [(dd["points_all"][0], "s", "Codex criterion"),
                         (dd["points_all"][1], "D", "TM-score \u2265 0.5")]:
        ax.plot(1 - pt["spec"], pt["sens"], mk_, ms=4.5, color=F.RED, mew=0, label="_" + lab)
        pts_h.append(Line2D([], [], marker=mk_, ls="", ms=4.5, color=F.RED, mew=0, label=lab))
    ax.plot([0, 1], [0, 1], lw=0.7, color=F.METHOD["other"], ls=(0, (4, 3)))
    ax.set_xlim(-0.02, 1.02); ax.set_ylim(-0.02, 1.02)
    ax.set_xlabel("1 \u2212 specificity")
    ax.set_ylabel("Sensitivity")
    # 圖例放左上（曲線上方的空白），說明文字放右下（對角線下方的空白），兩者都不壓曲線
    leg = ax.legend(frameon=False, loc="upper left", fontsize=6.0, handletextpad=0.3,
                    borderaxespad=0.2, title="AUC (95% CI)", title_fontsize=6.0, labelspacing=0.35,
                    alignment="left", handlelength=1.0)
    leg2 = ax.legend(handles=pts_h, frameon=False, loc="lower right", fontsize=6.0,
                     handletextpad=0.4, bbox_to_anchor=(1.0, 0.36), borderaxespad=0.2)
    ax.add_artist(leg)
    jt = rb[("joint model", "TM (rank, SD)")]
    js = rb[("joint model", "sequence identity (rank, SD)")]
    cd = json.load(open(OUT / "coef_diff.json"))
    lo_, hi_ = cd["diff_ci"]
    note = ax.annotate("No difference detected.\n"
                "Joint model, homologous pairs,\n"
                f"per SD of rank: TM-score {jt['estimate']:.2f},\n"
                f"identity {js['estimate']:.2f}; difference {cd['diff']:.2f}\n"
                f"(95% CI {lo_:.2f} to {hi_:.2f})".replace("-0", "\u22120").replace("-1", "\u22121"),
                xy=(0.98, 0.03), xycoords="axes fraction", ha="left", va="bottom",
                fontsize=6.0)
    # 右對齊整塊但每行靠左：先量寬度再移
    fig.canvas.draw()
    bb = note.get_window_extent().transformed(ax.transAxes.inverted())
    note.set_position((0.98 - bb.width, 0.03))
    # 程式檢查：圖例與說明文字都不可壓到任何一條曲線（含對角線）
    fig.canvas.draw()
    for art in [*leg.get_texts(), leg.get_title(), *leg2.get_texts(), note]:   # 逐行字，不用整個圖例外框
        box = art.get_window_extent()
        for ln in ax.lines:
            xy = ax.transData.transform(np.column_stack(ln.get_data()))
            xs = np.linspace(xy[:, 0].min(), xy[:, 0].max(), 400) if len(xy) == 2 else None
            if xs is not None:
                xy = np.column_stack([xs, np.interp(xs, xy[:, 0], xy[:, 1])])
            inside = ((xy[:, 0] > box.x0) & (xy[:, 0] < box.x1) &
                      (xy[:, 1] > box.y0) & (xy[:, 1] < box.y1))
            if inside.any() and not ln.get_label().startswith("_"):
                pts = ax.transData.inverted().transform(xy[inside]); raise SystemExit(f"Figure 3D：{type(art).__name__} 壓到曲線 {ln.get_label()} {pts[[0, -1]]} {box.transformed(ax.transAxes.inverted())}")
    F.bare(ax); F.panel(ax, "D", dx=-36)
    overlap_check(fig)
    data_overlap_check(fig, skip=(fig.axes[1],))   # 色條
    F.save(fig, "Figure_3")


# ============================================================ Figure 4（生物學）
def _overlay(ax, a, b, meta, labels):
    """兩個 AlphaFold 模型的 cartoon 疊合（overlay_render.py 用 PyMOL 畫好的 PNG）。

    疊合矩陣是 overlay.py 的 foldseek TM-align 結果（已驗 RMSD），PNG 只是把它畫成 cartoon。
    這裡把白邊裁掉再擺進格子。"""
    img = plt.imread(OUT / f"overlay_{a}_{b}.png")
    ink = img[..., :3].min(-1) < 0.97          # 白底（不透明背景）以外的像素
    ys, xs = np.where(ink)
    img = img[ys.min():ys.max() + 1, xs.min():xs.max() + 1]
    ax.imshow(img, interpolation="lanczos")
    ax.set_anchor("N")      # C、D 裁切後長寬比不同，靠上對齊才會讓兩個標題同高
    ax.set_axis_off()
    m = meta[f"{a}|{b}"]
    ax.set_title(f"{labels[0]} \u2013 {labels[1]}\n", fontsize=7.4, pad=2)
    ax.text(0.5, 1.0, f"TM-score {m['tm']:.2f}   identity {m['seqid80']:.0f}%   "
                      f"OR {m['odds_ratio']:.1f}", transform=ax.transAxes,
            ha="center", va="bottom", fontsize=6.6, color=F.GREY)


def figure4():
    import networkx as nx
    c = core()
    meta = json.load(open(OUT / "overlay.json"))

    fig = plt.figure(figsize=(F.WIDTH_IN, 6.0))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 0.62], hspace=0.30, wspace=0.30,
                          left=0.03, right=0.95, top=0.98, bottom=0.03)

    # --- A 網路：非同源的邊只留極淡的底色，彩色節點才是重點
    ax = fig.add_subplot(gs[0, 0])
    strong = c[c.z > 3]
    G = nx.Graph()
    for t in strong.itertuples():
        G.add_edge(t.a, t.b, hom=bool(t.homologous))
    fam = families.table()
    fam = fam[fam.status == "ok"].drop_duplicates("component").set_index("component").family
    # 以家族為種子的起始位置，spring 只做微調：否則同家族會被非同源的邊拉散
    order = list(FAM_COLORS) + ["other"]
    ang = {f: 2 * np.pi * i / len(order) for i, f in enumerate(order)}
    rng = np.random.default_rng(3)
    init = {n: np.array([np.cos(ang.get(fam.get(n), ang["other"])),
                         np.sin(ang.get(fam.get(n), ang["other"]))]) * 1.0
            + rng.normal(0, 0.08, 2) for n in G.nodes()}
    pos = nx.spring_layout(G, pos=init, seed=7, k=0.30, iterations=60)
    eh = [(u, v) for u, v, d in G.edges(data=True) if d["hom"]]
    en = [(u, v) for u, v, d in G.edges(data=True) if not d["hom"]]
    nx.draw_networkx_edges(G, pos, edgelist=en, ax=ax, edge_color="#e9edf0", width=0.3,
                           alpha=0.5)
    nx.draw_networkx_edges(G, pos, edgelist=eh, ax=ax, edge_color=F.NAVY, width=0.8,
                           alpha=0.75)
    cols = [FAM_COLORS.get(fam.get(n), OTHER) for n in G.nodes()]
    nx.draw_networkx_nodes(G, pos, ax=ax, node_size=22, node_color=cols, linewidths=0)
    ax.axis("off")
    handles = [Line2D([], [], marker="o", ls="", ms=4, color=v, label=k.split(" (")[0])
               for k, v in FAM_COLORS.items()]
    handles += [Line2D([], [], marker="o", ls="", ms=4, color=OTHER, label="other family"),
                Line2D([], [], color=F.NAVY, lw=1.1, label="homologous pair"),
                Line2D([], [], color="#dfe4e8", lw=1.1, label="non-homologous pair")]
    ax.legend(handles=handles, frameon=False, fontsize=5.8, loc="upper center",
              bbox_to_anchor=(0.5, 0.06), ncol=3, handletextpad=0.3, labelspacing=0.3,
              columnspacing=0.8, borderaxespad=0)
    n_nonhom = len(en)
    ax.annotate(f"{n_nonhom:,} of {len(G.edges()):,} edges join components that are not "
                f"homologous", xy=(0.5, 1.0), xycoords="axes fraction", ha="center", va="top",
                fontsize=6.2, color=F.GREY)
    F.panel(ax, "A", dx=-4, dy=-6)

    # --- B 家族 × 家族，格子裡直接寫數值（純色階在對角線會飽和成同一個紅）
    ax = fig.add_subplot(gs[0, 1])
    fams = list(FAM_COLORS)
    M = np.full((len(fams), len(fams)), np.nan)
    for i, fa in enumerate(fams):
        for j, fb in enumerate(fams):
            g = c[((c.fam_a == fa) & (c.fam_b == fb)) | ((c.fam_a == fb) & (c.fam_b == fa))]
            if len(g) >= 3:
                M[i, j] = g.beta.median()
    im = ax.imshow(M, cmap="RdYlBu_r", norm=TwoSlopeNorm(vmin=-2, vcenter=0, vmax=5))
    for i in range(len(fams)):
        for j in range(len(fams)):
            if not np.isnan(M[i, j]):
                v = M[i, j]
                ax.text(j, i, f"{v:.1f}", ha="center", va="center", fontsize=5.0,
                        color="white" if v > 3.2 or v < -1.4 else F.INK)
    short = [f.split(" (")[0] for f in fams]
    ax.set_xticks(range(len(fams)))
    ax.set_xticklabels(short, rotation=45, ha="right", rotation_mode="anchor", fontsize=5.8)
    ax.set_yticks(range(len(fams))); ax.set_yticklabels(short, fontsize=5.8)
    ax.tick_params(length=0)
    for sp in ax.spines.values():
        sp.set_visible(False)
    cb = fig.colorbar(im, ax=ax, pad=0.03, fraction=0.045, extend="max")   # 對角線有 > 5 的格
    cb.set_label("Median log odds ratio", fontsize=6.6)
    cb.ax.tick_params(labelsize=6.2)
    cb.outline.set_linewidth(0.5)
    F.panel(ax, "B", dx=-52, dy=-6)

    # --- C、D 兩組疊合：同摺疊但序列測不到 vs 同家族但摺疊已分化
    axc = fig.add_subplot(gs[1, 0])
    _overlay(axc, "Cry_j_1", "Pla_a_2", meta, ("Cry j 1", "Pla a 2"))
    F.panel(axc, "C", dx=4, dy=6)
    axd = fig.add_subplot(gs[1, 1])
    _overlay(axd, "Ara_h_2", "Api_g_2", meta, ("Ara h 2", "Api g 2"))
    F.panel(axd, "D", dx=4, dy=6)
    overlap_check(fig, skip=(ax, axc, axd, fig.axes[0]))
    # 上一版 B 的橫軸標籤與 A 的圖例被下排蓋掉／切掉而 overlap_check 沒抓到（那兩格被 skip）：
    # 這裡另外要求它們整個落在畫布內，且不碰到 C、D 的標題列
    rr = fig.canvas.get_renderer()
    W = fig.bbox
    lower = [t.get_window_extent(rr) for a_ in (axc, axd) for t in [a_.title, *a_.texts]]
    lower += [a_.get_window_extent(rr) for a_ in (axc, axd)]
    upper = [t for t in ax.get_xticklabels() if t.get_text()] + [fig.axes[0].get_legend()]
    for t in upper:
        bb = t.get_window_extent(rr)
        if bb.x0 < W.x0 - 1 or bb.x1 > W.x1 + 1 or bb.y0 < W.y0 - 1 or bb.y1 > W.y1 + 1:
            raise SystemExit(f"Figure 4 的 {getattr(t, 'get_text', lambda: 'legend')()!r} 超出畫布")
        if any(bb.overlaps(o) for o in lower):
            raise SystemExit(f"Figure 4 的 {getattr(t, 'get_text', lambda: 'legend')()!r} 碰到下排")
    F.save(fig, "Figure_4")


# ============================================================ 補充圖
def figure_s1():
    """敏感度分析森林圖（含節點自助法信賴區間）＋留一同源群。"""
    ci = json.load(open(OUT / "rho_ci.json"))
    loo = pd.read_csv(OUT / "robust_leave_one_group.csv")
    rows = [("Homologous pairs, main analysis", "homologous"),
            ("  ISAC v1 only", "isac_v1"),
            ("  ISAC v2 only", "isac_v2"),
            ("  ALEX2 only", "alex"),
            ("  Source species in different orders", "different_order"),
            ("  Both AlphaFold models pLDDT \u2265 70", "plddt70"),
            ("  Below the sequence criterion", "below_codex"),
            ("Non-homologous pairs, main analysis", "nonhomologous"),
            ("  Excluding fragment sequences", "nonhomologous_nofrag"),
            ("  Curated family labels instead of Pfam/SCOP", "nonhomologous_curated")]
    fig = plt.figure(figsize=(F.WIDTH_IN, 3.4))
    gs = fig.add_gridspec(1, 2, width_ratios=[2.05, 1.0], wspace=0.10,
                          left=0.33, right=0.95, bottom=0.17, top=0.96)

    ax = fig.add_subplot(gs[0, 0])
    y = np.arange(len(rows))[::-1]
    for yy, (lab, k) in zip(y, rows):
        v = ci[k]
        col = F.NAVY if not lab.startswith(" ") else F.METHOD["other"]
        ax.plot([v["lo"], v["hi"]], [yy, yy], lw=1.0, color=col, zorder=2)
        ax.scatter([v["rho"]], [yy], s=16, color=col, lw=0, zorder=3)
        ax.annotate(f"n = {v['n']:,}", xy=(1.28, yy), ha="right", va="center",
                    fontsize=6.6, color=F.GREY)
    ax.axvline(0, lw=0.8, color=F.METHOD["other"], ls=(0, (4, 3)))
    ax.set_yticks(y); ax.set_yticklabels([r[0] for r in rows], fontsize=7.0)
    ax.set_xlim(-0.25, 1.30); ax.set_ylim(-0.7, len(rows) - 0.3)
    ax.set_xticks([-0.2, 0, 0.2, 0.4, 0.6, 0.8, 1.0])
    ax.set_xlabel("Spearman \u03c1 between TM-score and\nco-sensitisation (95% CI)",
                  fontsize=8)
    F.bare(ax); F.panel(ax, "A", dx=-4, dy=6)

    ax = fig.add_subplot(gs[0, 1])
    l = loo.sort_values("rho_tm")
    ax.scatter(l.rho_tm, np.arange(len(l)), s=11, color=F.NAVY, lw=0)
    full = ci["homologous"]["rho"]
    ax.axvline(full, lw=0.9, color=F.RED, ls=(0, (3, 2)))
    ax.set_yticks([]); ax.set_ylim(-1.6, len(l) + 0.6)
    ax.set_xlim(0.6, 0.9); ax.set_xticks([0.6, 0.7, 0.8, 0.9])
    ax.set_xlabel("\u03c1 after leaving out one\nhomology group", fontsize=8)
    ax.annotate(f"{len(l)} groups\n{l.rho_tm.min():.2f} to {l.rho_tm.max():.2f}",
                xy=(0.04, 0.97), xycoords="axes fraction", va="top", fontsize=6.8)
    ax.annotate("all pairs", xy=(full, -1.3), xytext=(3, 0), textcoords="offset points",
                ha="left", va="bottom", fontsize=6.4, color=F.RED)
    F.bare(ax); F.panel(ax, "B", dx=-4, dy=6)
    overlap_check(fig)
    data_overlap_check(fig)
    F.save(fig, "Figure_S1")


def figure_s2():
    """AlphaFold 模型信心。"""
    pl = json.load(open(OUT / "plddt.json"))
    cp = pd.read_csv(OUT / "component_plddt.csv")
    c = core()
    lut = cp.set_index("component").plddt
    hom = c[c.homologous].copy()
    hom["plddt_min"] = np.minimum(hom.a.map(lut), hom.b.map(lut))
    fig, axes = plt.subplots(1, 2, figsize=(F.WIDTH_IN, 2.5))
    ax = axes[0]
    lo_ = 2.5 * np.floor(cp.plddt.min() / 2.5)
    hh_, _, _ = ax.hist(cp.plddt, bins=np.arange(lo_, 102.5, 2.5), color=F.NAVY, alpha=0.85)
    assert hh_.sum() == len(cp), "直方圖漏畫成分"
    ax.axvline(70, lw=0.9, color=F.RED, ls=(0, (3, 2)))
    ax.set_xlabel("Mean pLDDT of the AlphaFold model")
    ax.set_ylabel("Components")
    ax.annotate(f"{pl['n_components_below']} of {pl['n_components']} below 70",
                xy=(0.03, 0.95), xycoords="axes fraction", va="top", fontsize=7)
    F.bare(ax); F.panel(ax, "A", dx=-32)
    ax = axes[1]
    ax.scatter(hom.plddt_min, hom.tm, s=7, lw=0, color=F.METHOD["main"], alpha=0.5, zorder=3)
    ax.set_ylim(0.2, 1.42)
    ax.set_yticks([0.2, 0.4, 0.6, 0.8, 1.0])
    ax.plot([70, 70], [0.2, 1.03], lw=0.9, color=F.RED, ls=(0, (3, 2)), zorder=1)
    ax.set_xlabel("Lower pLDDT of the two models")
    ax.set_ylabel("TM-score")
    ax.annotate(f"ρ = {pl['rho_plddt_tm']:.2f}\nρ(TM, co-sensitisation) "
                f"{pl['rho_all']:.2f} → {pl['rho_high_confidence']:.2f}\nwhen both ≥ 70",
                xy=(0.03, 0.95), xycoords="axes fraction", va="top", fontsize=7)
    F.bare(ax); F.panel(ax, "B", dx=-32)
    fig.tight_layout(w_pad=2.2)
    overlap_check(fig)
    data_overlap_check(fig)
    F.save(fig, "Figure_S2")


def figure_s3():
    """排列檢定的虛無分布。"""
    b = json.load(open(OUT / "bands.json"))
    nh = np.load(OUT / "null_homologous.npy")
    nn = np.load(OUT / "null_nonhomologous.npy")
    fig, axes = plt.subplots(1, 2, figsize=(F.WIDTH_IN, 2.4))
    for ax, null, obs, lab, letter in [
            (axes[0], nh, b["obs_homologous"], "Homologous pairs", "A"),
            (axes[1], nn, b["obs_nonhomologous"], "Non-homologous pairs", "B")]:
        hh, _, _ = ax.hist(null, bins=40, color=OTHER, lw=0)
        top = hh.max() * 1.5                         # 上方留白給註記，字不壓直方圖
        ax.set_ylim(0, top)
        ax.axvline(obs, lw=1.3, color=F.RED)
        ax.annotate(f"observed\nρ = {obs:.3f}", xy=(obs, top * 0.98),
                    ha="right" if obs > np.median(null) else "left", va="top", fontsize=7,
                    color=F.RED, xytext=(-4 if obs > np.median(null) else 4, 0),
                    textcoords="offset points")
        ax.set_xlabel("Spearman ρ under permuted allergen labels")
        ax.set_ylabel("Permutations")
        ax.set_title(lab, fontsize=7.6, pad=3)
        F.bare(ax); F.panel(ax, letter, dx=-34)
    fig.tight_layout(w_pad=2.2)
    overlap_check(fig)
    data_overlap_check(fig)
    F.save(fig, "Figure_S3")


def figure_s4():
    """三個平台各自估計的一致性。"""
    bc = pd.read_csv(OUT / "cosens_by_chip.csv")
    bc = bc[bc.method == "mle"]
    c = core()[["a", "b", "tm", "homologous"]]
    fig, axes = plt.subplots(1, 3, figsize=(F.WIDTH_IN, 2.5))
    for ax, chip, letter in zip(axes, F.CHIP_ORDER, "ABC"):
        m = c.merge(bc[bc.chip == chip][["a", "b", "beta"]], on=["a", "b"])
        h = m[m.homologous]
        ax.scatter(m.tm[~m.homologous], m.beta[~m.homologous], s=3, lw=0, alpha=0.2,
                   color=F.METHOD["other"], rasterized=True)
        ax.scatter(h.tm, h.beta, s=7, lw=0, alpha=0.75, color=F.CHIP[chip])
        ax.axhline(0, lw=0.7, color=F.METHOD["other"], ls=(0, (4, 3)))
        ax.set_xlim(0.15, 1.0); ax.set_ylim(-6.5, 10)
        ax.set_xlabel("TM-score")
        if letter == "A":
            ax.set_ylabel("Adjusted log odds ratio")
        from scipy import stats as st
        r = st.spearmanr(h.tm, h.beta).statistic
        ax.set_title(f"{F.CHIP_LABEL[chip]}  (n = {len(h)}, ρ = {r:.2f})", fontsize=7.4, pad=3)
        F.bare(ax); F.panel(ax, letter, dx=-30)
    fig.tight_layout(w_pad=1.6)
    overlap_check(fig)
    data_overlap_check(fig)
    F.save(fig, "Figure_S4")


if __name__ == "__main__":
    import sys
    todo = sys.argv[1:] or ["1", "2", "3", "4", "s1", "s2", "s3", "s4"]
    fn = {"1": figure1, "2": figure2, "3": figure3, "4": figure4,
          "s1": figure_s1, "s2": figure_s2, "s3": figure_s3, "s4": figure_s4}
    for k in todo:
        print(k, flush=True)
        fn[k]()
