#!/usr/bin/env python3
"""绘制同一初态下原、修正FG级数相对二体加J2数值预报的位置差。"""
import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--csv", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path, help="输出PNG、PDF、SVG的文件前缀。")
    parser.add_argument("--qa-dir", type=Path, help="可选布局检查目录；需要nature-figure的审查脚本在PYTHONPATH中。")
    args = parser.parse_args()
    data = np.loadtxt(args.csv, delimiter=",", skiprows=1)
    assert data.shape == (601, 4) and np.isfinite(data).all()
    assert np.all(np.diff(data[:, 0]) > 0) and np.all(data[:, 1:] >= 0)

    # 共同初态对应零误差，纵轴在微米以下使用线性区间，从而保留历元处的数据。
    shown = data[data[:, 0] >= 0]
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial Unicode MS", "Hiragino Sans GB", "DejaVu Sans"],
        "font.size": 9, "axes.labelsize": 9, "xtick.labelsize": 8, "ytick.labelsize": 8,
        "legend.fontsize": 9, "legend.frameon": False, "pdf.fonttype": 42,
        "svg.fonttype": "none", "axes.spines.top": False, "axes.spines.right": False,
        "axes.linewidth": 0.7, "axes.unicode_minus": False, "figure.facecolor": "white",
    })
    fig, ax = plt.subplots(figsize=(7.2, 4.3))
    fig.subplots_adjust(left=0.13, right=0.975, bottom=0.16, top=0.73)
    ax.plot(shown[:, 0], shown[:, 1], color="#A45A45", linestyle="--", linewidth=1.7, label="原 FG 级数")
    ax.plot(shown[:, 0], shown[:, 2], color="#356CA5", linewidth=1.7, label="修正 FG 级数")
    ax.set(xlim=(0, 315), ylim=(0, 1000), xlabel="距展开历元的时间 / s",
           ylabel="相对数值预报的三维位置差 / m")
    ax.set_yscale("symlog", linthresh=1e-6, linscale=0.5)
    ax.set_yticks([0, 1e-6, 1e-4, 1e-2, 1, 100])
    ax.set_yticklabels(["0", "10⁻⁶", "10⁻⁴", "10⁻²", "1", "100"])
    ax.grid(axis="y", color="#E7E9EC", linewidth=0.6, zorder=0)
    ax.legend(loc="lower left", bbox_to_anchor=(0, 1.02), ncol=2, borderaxespad=0)
    fig.text(0.13, 0.94, "修正 Gz 的 J₂ 三次项后，预报更接近数值积分", fontsize=12, va="top")
    fig.text(0.13, 0.855, "a = 7000 km，e = 0，i = 55°；同一初态、相同二体 + J₂ 模型", fontsize=9, va="top")

    # 末点直接标注两条曲线的实际数值；完整负向时间序列仍保存在CSV中。
    for column, color in ((1, "#A45A45"), (2, "#356CA5")):
        value = shown[-1, column]
        ax.plot(300, value, "o", color=color, markersize=3.5)
        ax.annotate(f"{value:.2f} m", xy=(300, value), xytext=(-6, 9),
                    textcoords="offset points", ha="right", color=color, fontsize=9)

    if args.qa_dir:
        from audit_panel_alignment import require_matplotlib_panel_alignment
        args.qa_dir.mkdir(parents=True, exist_ok=True)
        require_matplotlib_panel_alignment(
            fig, json_out=args.qa_dir / "iod_j2_cubic_comparison.alignment.json",
            overlay_svg=args.qa_dir / "iod_j2_cubic_comparison.alignment.svg", strict=True)

    args.output.parent.mkdir(parents=True, exist_ok=True)

    fig.savefig(args.output.with_suffix(".png"), dpi=600)
    fig.savefig(args.output.with_suffix(".pdf"))
    fig.savefig(args.output.with_suffix(".svg"))
    svg = args.output.with_suffix(".svg")
    svg.write_text("\n".join(line.rstrip() for line in svg.read_text().splitlines()) + "\n")
    plt.close(fig)


if __name__ == "__main__":
    main()
