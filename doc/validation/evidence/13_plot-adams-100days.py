"""从双方完整 Cartesian 星历绘制百天位置差；输入位置单位为 km。"""
from pathlib import Path
import argparse
import hashlib
import json
import os
os.environ.setdefault('MPLCONFIGDIR', '/tmp/aster-matplotlib')
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--data-dir', type=Path, default=Path(__file__).resolve().parent)
parser.add_argument('--output-dir', type=Path)
args = parser.parse_args()
source = args.data_dir.resolve()
out = args.output_dir or source / 'figures'
out.mkdir(parents=True, exist_ok=True)

# 使用中文字体；曲线保留每一个规则网格采样点，最终时刻取各自补齐后的输出。
font_path = Path('/System/Library/Fonts/Supplemental/Arial Unicode.ttf')
if font_path.exists():
    font_manager.fontManager.addfont(str(font_path))
    plt.rcParams['font.family'] = font_manager.FontProperties(fname=str(font_path)).get_name()
plt.rcParams.update({'font.size': 10, 'axes.spines.top': False, 'axes.spines.right': False,
                     'pdf.fonttype': 42, 'svg.fonttype': 'none', 'axes.unicode_minus': False, 'path.simplify': False})
fig, ax = plt.subplots(figsize=(9.2, 5.8))
fig.subplots_adjust(left=0.12, right=0.96, bottom=0.15, top=0.78)
metadata = {'comparison': 'ASTER vs mentor NAO, Cartesian position norm in metres',
            'alignment': 'same nominal elapsed time; final sample uses final output of each program',
            'duration_days': 100, 'integrator': 'Adams PECE, 9 steps', 'curves': []}

for step, suffix, color in [(60, '', '#B86736'), (30, '_30s', '#3479AA'), (10, '_10s', '#25816D')]:
    paths = [source / ('aster' + suffix) / 'cart.txt', source / ('nao' + suffix) / 'cart.txt']
    a = np.loadtxt(paths[0]); b = np.loadtxt(paths[1])[:, 1:]
    n = 8640000 // step + 1
    assert len(a) in (n, n + 1) and len(b) in (n, n + 1)
    assert np.isfinite(a).all() and np.isfinite(b).all()
    assert np.all(np.diff(a[:, 0]) > 0) and np.all(np.diff(b[:, 0]) > 0)
    counts = [len(a), len(b)]
    # 补齐步与最后一个规则点只差微秒；用实际终点替换该规则点。
    endpoint_a, endpoint_b = a[-1].copy(), b[-1].copy()
    a, b = a[:n].copy(), b[:n].copy()
    a[-1], b[-1] = endpoint_a, endpoint_b
    delta = np.linalg.norm(a[:, 1:4] - b[:, 1:4], axis=1) * 1000
    day = np.arange(n) * step / 86400
    ax.plot(day, delta, color=color, linewidth=1.1, label=f'{step} s', rasterized=True)
    metadata['curves'].append({'step_s': step, 'source_rows': counts, 'plotted_points': n,
        'max_m': float(delta.max()), 'max_day': float(day[delta.argmax()]), 'final_m': float(delta[-1]),
        'max_printed_epoch_difference_s': float(np.max(np.abs(a[:, 0] - b[:, 0])) * 86400),
        'files': [{'path': str(p), 'sha256': hashlib.sha256(p.read_bytes()).hexdigest()} for p in paths]})

ax.set(xlim=(0, 100), ylim=(0, 2400), xlabel='预报时长 / 天', ylabel='ASTER 与导师程序的三维位置差 / m')
ax.set_xticks(np.arange(0, 101, 10))
ax.set_yticks(np.arange(0, 2401, 400))
ax.grid(axis='y', color='#DDE2E6', linewidth=0.6)
ax.set_axisbelow(True)
fig.text(0.12, 0.94, '100 天轨道预报：ASTER 与导师程序对照', fontsize=15, weight='bold')
fig.text(0.12, 0.885, '同一卫星与初值 · ADAMS PECE（9 步）· 各步长分别比较双方完整星历', fontsize=10, color='#424B55')
ax.legend(title='固定步长', loc='upper left', frameon=False, labelspacing=0.8)
fig.text(0.12, 0.035, '地球 80 阶、日月质点、阻力、光压及地影开启；ASTER 另含月球 20 阶非球形。', fontsize=9, color='#424B55')

stem = out / 'adams_100days_comparison'
fig.savefig(stem.with_suffix('.png'), dpi=300)
fig.savefig(stem.with_suffix('.pdf'), dpi=450)
fig.savefig(stem.with_suffix('.svg'), dpi=450)
plt.close(fig)
metadata['rendering'] = {'all_regular_samples_used': True, 'smoothing': False,
    'pdf': 'vector text and axes; all curve samples rasterized at 450 dpi', 'alignment_audit': 'not applicable: one panel'}
(out / 'plot_metadata.json').write_text(json.dumps(metadata, ensure_ascii=False, indent=2))
print(json.dumps(metadata['curves'], ensure_ascii=False, indent=2))
