#!/usr/bin/env python3
"""绘制四套程序或潮汐模型相对当前 ASTER 的三维位置差。"""
import argparse
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--csv', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True, help='输出文件前缀，保存 PNG、PDF、SVG。')
    parser.add_argument('--qa-dir', type=Path, help='保存单图布局检查记录。')
    args = parser.parse_args()

    # CSV 保留所有历元，时间为初始事件后的小时数，位置差单位为 m。
    data = np.loadtxt(args.csv, delimiter=',', skiprows=1)
    assert data.shape == (8641, 5) and np.isfinite(data).all()
    assert np.all(np.diff(data[:, 0]) > 0) and np.all(data[:, 1:] >= 0)
    plt.rcParams.update({
        'font.family': 'sans-serif',
        'font.sans-serif': ['Arial Unicode MS', 'Hiragino Sans GB', 'DejaVu Sans'],
        'font.size': 9,
        'axes.spines.top': False,
        'axes.spines.right': False,
        'axes.linewidth': 0.7,
        'axes.unicode_minus': False,
        'pdf.fonttype': 42,
        'svg.fonttype': 'none',
        'legend.frameon': False,
        'figure.facecolor': 'white',
    })

    # 单张主图同时显示米级和微米级差异，零点附近使用线性刻度保留共同初态。
    fig = plt.figure(figsize=(7.2, 4.5))
    ax = fig.add_axes([0.12, 0.16, 0.85, 0.57])
    ax.set_yscale('symlog', linthresh=1e-6, linscale=0.5, base=10)
    specifications = [
        ('NAO 原模型（重复永久潮）', '#B7473E', '-', 1.5),
        ('NAO 扣除重复永久潮', '#356CA5', '--', 1.5),
        ('Orekit 二阶对齐模型', '#278A70', '-', 1.1),
        ('Orekit IERS 2010', '#BC812D', '-.', 1.5),
    ]

    for column, (label, color, linestyle, width) in enumerate(specifications, start=1):
        ax.plot(data[:, 0], data[:, column], label=label, color=color,
                linestyle=linestyle, linewidth=width)

    # 采用同一组时刻和单位，图例逐项说明程序或模型的实际设置。
    ax.set_xlim(0, 24)
    ax.set_ylim(-8e-8, 20)
    ax.set_xticks(np.arange(0, 25, 4))
    ax.set_yticks([0, 1e-6, 1e-4, 1e-2, 1, 10])
    ax.set_yticklabels(['0', r'$10^{-6}$', r'$10^{-4}$', r'$10^{-2}$', '1', '10'])
    ax.set_xlabel('预报时间 / h', labelpad=7)
    ax.set_ylabel('相对 ASTER 的位置差 / m', labelpad=8)
    ax.tick_params(axis='both', which='major', labelsize=9, length=3, width=0.7)
    ax.tick_params(axis='y', which='minor', left=False)
    ax.grid(axis='y', which='major', color='#E0E4E8', linewidth=0.5)
    ax.set_axisbelow(True)

    fig.text(0.12, 0.965, '地球固体潮轨道预报对照', fontsize=12, weight='bold', va='top')
    fig.text(0.12, 0.907, 'GGM03C 20×20 + 固体潮 · 24 h 预报 · 基线：当前 ASTER', fontsize=8.5, color='#555F68')
    fig.legend(*ax.get_legend_handles_labels(), loc='upper left', bbox_to_anchor=(0.11, 0.877),
               ncol=2, fontsize=8.5, columnspacing=2.0, handlelength=2.7, labelspacing=0.8)
    fig.text(0.12, 0.025, '纵轴在 1 μm 以内为线性刻度，其余为对数刻度。', fontsize=7.5, color='#59636D')

    # 单图布局记录为不适用多面板对齐；以最终导出尺寸检查文字和图例位置。
    fig.canvas.draw()

    if args.qa_dir:
        from audit_panel_alignment import require_matplotlib_panel_alignment
        args.qa_dir.mkdir(parents=True, exist_ok=True)
        require_matplotlib_panel_alignment(fig, json_out=args.qa_dir / 'alignment.json',
            overlay_svg=args.qa_dir / 'alignment.svg', tolerance_pt=1.5, strict=True)

    args.output.parent.mkdir(parents=True, exist_ok=True)

    fig.savefig(f'{args.output}.png', dpi=600, facecolor='white')
    fig.savefig(f'{args.output}.pdf', facecolor='white')
    fig.savefig(f'{args.output}.svg', facecolor='white')

    plt.close(fig)


if __name__ == '__main__':
    main()
