#!/usr/bin/env python3
"""检查双向单程测距样例、噪声设置和 Tracer 错误出口。"""
import argparse
import csv
import json
import math
from pathlib import Path
import subprocess
import tempfile

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--binary', type=Path, required=True)
parser.add_argument('--case', type=Path, required=True)
parser.add_argument('--output', type=Path, required=True)
args = parser.parse_args()
args.output.mkdir(parents=True, exist_ok=True)
args.output = Path(tempfile.mkdtemp(prefix="run-", dir=args.output))
base = (args.case / 'aster.settings').read_text()

# 每种设置使用独立输出目录，保留运行日志供定位。
for name, changes, expected_success in [
    ('basic', {}, True),
    ('bias', {'systematic_bias_m': '2.5'}, True),
    ('noise_a', {'noise_sigma_m': '1.0', 'noise_seed': '123'}, True),
    ('noise_b', {'noise_sigma_m': '1.0', 'noise_seed': '123'}, True),
    ('no_gravity', {'gravity_delay': '0'}, True),
    ('missing_input', {'file_user_state': 'missing_state.txt'}, False),
    ('outside_epoch', {'file_receive_events': str((args.output / 'outside.txt').resolve())}, False),
]:
    settings = dict(line.split('=', 1) for line in base.splitlines() if '=' in line and not line.startswith('#'))
    settings = {key.strip(): value.strip() for key, value in settings.items()}
    settings.update(changes)
    settings['directory_output'] = str((args.output / name).resolve())
    for key in ['file_navigation_state', 'file_user_state', 'file_receive_events', 'file_JPLEPH']:
        settings[key] = str((args.case / settings[key]).resolve())
    (args.output / 'outside.txt').write_text('2029 5 20 5 59 0 12 901\n')
    config = args.output / f'{name}.settings'
    config.write_text('\n'.join(f'{k} = {v}' for k, v in settings.items()) + '\n')
    run = subprocess.run([str(args.binary.resolve()), '-I', str(config.resolve())], capture_output=True, text=True)
    (args.output / f'{name}.stderr').write_text(run.stderr)
    assert (run.returncode == 0) == expected_success, (name, run.stderr)
    if not expected_success:
        assert 'main' in run.stderr and '[ERROR]' in run.stderr, run.stderr
        if name == 'outside_epoch':
            assert 'InterpolateState' in run.stderr
        continue
    rows = list(csv.DictReader(line for line in (args.output / name / 'components.csv').read_text().splitlines() if not line.startswith('#')))
    assert len(rows) == 2
    for row in rows:
        assert all(math.isfinite(float(v)) for v in row.values())
        assert abs(float(row['light_time_s']) - float(row['geometric_m']) / 299792458.0 - float(row['gravity_s'])) <= 1.1e-12
        assert abs(float(row['clock_m']) - 299792458.0 * (float(row['rx_clock_s']) - float(row['tx_clock_s']))) < 1e-7
        assert abs(float(row['pseudorange_m']) - float(row['theoretical_pseudorange_m']) - float(row['systematic_bias_m']) - float(row['noise_m'])) < 1e-7
        if name == 'no_gravity':
            assert float(row['gravity_s']) == 0
        if name == 'bias':
            assert float(row['systematic_bias_m']) == 2.5
    if name == 'basic':
        assert [(row['tx'], row['rx']) for row in rows] == [('12', '901'), ('901', '12')]
        reference = list(csv.DictReader((args.case / 'reference.csv').open()))
        for row, expected in zip(rows, reference):
            for key in ['geometric_m', 'pseudorange_m', 'clock_m']:
                assert abs(float(row[key]) - float(expected[key])) < 1e-4, (key, row, expected)
assert (args.output / 'noise_a/observations.txt').read_bytes() == (args.output / 'noise_b/observations.txt').read_bytes()
(args.output / 'verification.json').write_text(json.dumps({'passed': True, 'cases': 7}, indent=2) + '\n')
print('Passed: bidirectional sample, bias, reproducible noise, gravity switch and Tracer failures.')
