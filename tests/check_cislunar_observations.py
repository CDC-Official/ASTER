#!/usr/bin/env python3
"""核对 448 条 LHH 算例，保存分方向统计、前两条样例及逐条差异。"""
import argparse
import csv
from datetime import datetime, timedelta
import hashlib
import json
import math
from pathlib import Path

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--aster', type=Path, required=True)
parser.add_argument('--reference', type=Path, required=True)
parser.add_argument('--diagnostic', type=Path, required=True)
parser.add_argument('--diagnostic-observations', type=Path, required=True)
parser.add_argument('--output', type=Path, required=True)
args = parser.parse_args()

# 诊断版只增加内部量输出，先确认原格式结果逐字节不变。
assert args.reference.read_bytes() == args.diagnostic_observations.read_bytes()
lines = args.aster.read_text().splitlines()
metadata = dict(line[2:].split('=', 1) for line in lines if line.startswith('# '))
rows = list(csv.DictReader(line for line in lines if not line.startswith('#')))
reference = [line.split() for line in args.reference.read_text().splitlines() if line.strip()]
diagnostics = [[float(x) for x in line.split()] for line in args.diagnostic.read_text().splitlines()]
assert len(rows) == len(reference) == len(diagnostics) == 448
assert metadata['jpl_de'] == '440'
origin = datetime(1858, 11, 17) + timedelta(days=float(metadata['reference_mjd_tdb_integer'])
                                        + float(metadata['reference_mjd_tdb_fraction']))

for row, old in zip(rows, reference):
    assert [row['type'], row['tx'], row['tx_antenna'], row['rx'], row['rx_antenna']] == [old[0], *old[2:6]]
    tag_s = (datetime.fromisoformat(old[1]) - origin).total_seconds()
    assert abs(float(row['clock_tag_s']) - tag_s) < 1e-9
    assert all(math.isfinite(float(value)) for key, value in row.items() if key != 'type')
    assert float(row['tx_tdb_s']) < float(row['rx_tdb_s'])
    assert abs(float(row['pseudorange_m']) - float(row['propagation_m']) - float(row['clock_m'])) < 1e-7
    assert abs(float(row['rx_tdb_s']) - float(row['tx_tdb_s']) - float(row['propagation_m']) / 299792458.0) < 3e-11
    assert int(row['rx_iterations']) <= 10 and int(row['tx_iterations']) <= 10
    assert float(row['rx_update_s']) <= 1e-12 and float(row['tx_update_s']) <= 1e-12

# LHH 内部按接收钟时稳定排序，原格式输出再恢复时隙表顺序。
ordered = sorted(rows, key=lambda row: float(row['clock_tag_s']))
fields = {'tx_tdb_s': (1, 86400.0), 'rx_tdb_s': (2, 86400.0), 'geometric_m': (3, 1.0),
          'gravity_s': (4, 1.0), 'propagation_m': (5, 1.0), 'gravity_m': (6, 1.0),
          'tx_physical_s': (7, 1.0), 'tx_relativistic_s': (8, 1.0),
          'rx_physical_s': (9, 1.0), 'rx_relativistic_s': (10, 1.0), 'pseudorange_m': (11, 1.0)}
differences = []
samples = []
for row, diagnostic in zip(ordered, diagnostics):
    assert abs(float(row['clock_tag_s']) - diagnostic[0] * 86400.0) < 1e-9
    difference = {'record': int(row['record']), 'tx': int(row['tx']), 'rx': int(row['rx'])}
    lhh = {name: diagnostic[index] * scale for name, (index, scale) in fields.items()}
    for name in fields:
        difference[name] = float(row[name]) - lhh[name]
    differences.append(difference)
    if int(row['record']) in (1, 2):
        samples.append({'record': int(row['record']), 'aster': row, 'lhh': lhh, 'aster_minus_lhh': difference})

differences.sort(key=lambda row: row['record'])
samples.sort(key=lambda row: row['record'])
summary = {}
for direction in ('all', 'BDS_to_901', '901_to_BDS'):
    selected = [row for row in differences if direction == 'all'
                or (direction == 'BDS_to_901' and row['rx'] == 901)
                or (direction == '901_to_BDS' and row['tx'] == 901)]
    summary[direction] = {'count': len(selected), 'metrics': {}}
    for name in fields:
        values = [row[name] for row in selected]
        worst = max(selected, key=lambda row: abs(row[name]))
        summary[direction]['metrics'][name] = {
            'max_abs': abs(worst[name]), 'rms': math.sqrt(sum(value * value for value in values) / len(values)),
            'worst_record': worst['record']}
assert summary['BDS_to_901']['count'] == summary['901_to_BDS']['count'] == 224
for name, metric in summary['all']['metrics'].items():
    tolerance = 1e-4 if name in ('geometric_m', 'propagation_m', 'pseudorange_m') else (
        1e-9 if name == 'gravity_m' else (1e-11 if name in ('tx_tdb_s', 'rx_tdb_s') else 1e-16))
    assert metric['max_abs'] <= tolerance, (name, metric, tolerance)

args.output.mkdir(parents=True, exist_ok=True)
with (args.output / 'differences.csv').open('w') as stream:
    writer = csv.DictWriter(stream, fieldnames=differences[0].keys())
    writer.writeheader()
    writer.writerows(differences)
report = {'metadata': metadata, 'diagnostic_preserves_original_output': True, 'checks_passed': True,
          'summary': summary, 'samples': samples,
          'sha256': {str(path): hashlib.sha256(path.read_bytes()).hexdigest()
                     for path in (args.aster, args.reference, args.diagnostic, args.diagnostic_observations)}}
(args.output / 'comparison.json').write_text(json.dumps(report, ensure_ascii=False, indent=2) + '\n')
print(json.dumps(summary['all'], ensure_ascii=False, indent=2))
