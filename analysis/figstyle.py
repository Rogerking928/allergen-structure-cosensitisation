"""本專案共用的繪圖設定。

沿用 Bcc 生合成島圖譜／嚴重氣喘那篇的版式：DejaVu Sans、細軸線、**無格線**，
以印出來的寬度作畫（180 mm = 7.09 in），存 700 dpi PNG。
**只存 PNG，不存 PDF**（使用者明確要求過兩次）。

過敏原代號（Bet v 1）照 WHO/IUIS 慣例用正體；只有物種學名要斜體。
"""
import os
import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
OUT = os.path.join(ROOT, "out")

WIDTH_IN = 7.09          # 180 mm
DPI = 700

NAVY = "#31485d"
RED = "#a4303f"
AMBER = "#b8791f"
GREEN = "#4f7d5b"
PURPLE = "#6b5b95"
TEAL = "#2f7d8c"
BROWN = "#7a5c3e"
PLUM = "#9c6b8a"
GREY = "#8a939c"
FAINT = "#c9d1d8"
INK = "#1f2933"

# 晶片平台：整份文件通用，讀者在第一張圖學到的對應之後不能改
CHIP = {"ISAC_V1": NAVY, "ISAC_V2": TEAL, "ALEX": AMBER}
CHIP_HATCH = {"ISAC_V1": "", "ISAC_V2": "///", "ALEX": "..."}
CHIP_LABEL = {"ISAC_V1": "ISAC v1", "ISAC_V2": "ISAC v2", "ALEX": "ALEX2"}
CHIP_ORDER = ["ISAC_V1", "ISAC_V2", "ALEX"]

# 方法對照（校正前 vs 校正後）用灰階＋形狀，不用顏色
METHOD = {"main": "#2b333b", "other": "#9aa3ab", "faint": "#c9d1d8"}

plt.rcParams.update({
    "font.family": ["DejaVu Sans"],
    "font.size": 9,
    "axes.titlesize": 9.5,
    "axes.labelsize": 9,
    "legend.fontsize": 8.5,
    "xtick.labelsize": 8.5,
    "ytick.labelsize": 8.5,
    "axes.unicode_minus": False,
    "axes.edgecolor": "#5b6b7a",
    "axes.labelcolor": INK,
    "text.color": INK,
    "xtick.color": "#5b6b7a",
    "ytick.color": "#5b6b7a",
    "figure.facecolor": "white",
    "hatch.linewidth": 0.6,
})


def bare(ax, keep=("bottom", "left")):
    for side, sp in ax.spines.items():
        sp.set_visible(side in keep)
    ax.grid(False)
    return ax


def panel(ax, letter, dx=-34, dy=8):
    ax.annotate(letter, xy=(0, 1), xycoords="axes fraction",
                xytext=(dx, dy), textcoords="offset points",
                fontsize=11, fontweight="bold", va="bottom", ha="left")


def check(fig, stem):
    """畫完後的自動檢查：撐寬是無聲失敗，必須擋下來（字會整頁變小）。"""
    probs = []
    fig.canvas.draw()
    bb = fig.get_tightbbox(fig.canvas.get_renderer())
    if bb.width > WIDTH_IN + 0.02:
        probs.append(f"tight bbox {bb.width:.2f} in > {WIDTH_IN} in")
    return probs


def save(fig, stem):
    os.makedirs(OUT, exist_ok=True)
    probs = check(fig, stem)
    if probs:
        raise SystemExit(f"{stem}: {probs[0]}")
    png = os.path.join(OUT, f"{stem}.png")
    fig.savefig(png, dpi=DPI, bbox_inches="tight", pad_inches=0.02, metadata={"Software": None})
    plt.close(fig)
    print(f"  fig -> {os.path.basename(png)}", flush=True)
