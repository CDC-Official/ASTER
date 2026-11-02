#!/usr/bin/env python3
"""运行七个既有非 IOD 算例，按现有跨平台物理量阈值比较科学输出。"""
import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import subprocess
import sys
import time

sys.path.insert(0, str(Path(__file__).resolve().parent))
from run_case import prepare

# 复用 doc/validation/evidence/03_linux-validation.json 的阈值。历元、记录编号和文本必须相同。
LIMITS = {'position_m': 1e-3, 'velocity_m_s': 1e-6, 'a_m': 1e-3,
          'eccentricity': 1e-10, 'inclination_parameters': 1e-10,
          'angle_deg': 1e-7, 'angle_residual_rad': 1e-10}
CASES = ['odn_batch', 'opn_cartesian', 'opn_nonsingular2', 'opn_cislunar',
         'convert_osculating_to_mean', 'convert_keplerian_to_cartesian',
         'convert_cartesian_to_keplerian']


def compare(actual, reference, case, name):
    a, b = actual.read_text().splitlines(), reference.read_text().splitlines()
    failures, maxima = [], {}
    if len(a) != len(b):
        return {'passed': False, 'reason': 'row count differs'}
    for row, (left, right) in enumerate(zip(a, b), 1):
        x, y = left.split(), right.split()
        if len(x) != len(y):
            failures.append([row, 'column count differs'])
            continue
        for col, (u, v) in enumerate(zip(x, y)):
            kind, scale = 'exact', 1.0
            if case == 'odn_batch':
                if col in (4, 5):
                    kind = 'angle_residual_rad'
            elif 'cartesian' in name or case == 'convert_keplerian_to_cartesian':
                if col > 0:
                    kind = 'position_m' if col <= 3 else 'velocity_m_s'
                    scale = 1.0 if case == 'convert_keplerian_to_cartesian' else 1000.0
            elif col == 1:
                kind, scale = 'a_m', 1000.0
            elif 'nsg2' in name:
                if col in (2, 3):
                    kind = 'eccentricity'
                elif col in (4, 5):
                    kind = 'inclination_parameters'
                elif col == 6:
                    kind = 'angle_deg'
            elif col == 2:
                kind = 'eccentricity'
            elif col > 2:
                kind = 'angle_deg'
            try:
                p, q = float(u), float(v)
            except ValueError:
                if u != v:
                    failures.append([row, col, u, v])
                continue
            if not math.isfinite(p) or not math.isfinite(q):
                if not (math.isnan(p) and math.isnan(q)) and p != q:
                    failures.append([row, col, u, v])
                continue
            error = abs(p - q) * scale
            if kind == 'angle_deg':
                error = abs(math.remainder(p - q, 360.0))
            maxima[kind] = max(maxima.get(kind, 0.0), error)
            if error > LIMITS.get(kind, 0.0):
                failures.append([row, col, error, kind])
    return {'passed': not failures, 'byte_equal': actual.read_bytes() == reference.read_bytes(),
            'maximum_errors': maxima, 'failures': failures[:10]}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--bin-dir', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--native-summary', type=Path)
    args = parser.parse_args()
    args.output = args.output.resolve()
    args.output.mkdir(parents=True, exist_ok=True)
    native = json.loads(args.native_summary.read_text()) if args.native_summary else None
    results = {'limits': LIMITS, 'cases': {}}
    passed = True
    for case in CASES:
        info, fixture, run, command = prepare(case, args.output / case)
        start = time.monotonic()
        process = subprocess.run([str(args.bin_dir.resolve() / info['target'])] + command,
                                 cwd=run, env=dict(os.environ, OMP_NUM_THREADS='2'),
                                 capture_output=True, timeout=120)
        elapsed = time.monotonic() - start
        (run / 'console.txt').write_bytes(process.stdout + process.stderr)
        result = {'directory': str(run), 'exit': process.returncode,
                  'elapsed_seconds': elapsed, 'files': {}}
        passed &= process.returncode == 0
        for name in info['references']:
            file = run / name
            reference = Path(native[case]['runs']['baseline']['run']) / name if native else fixture / 'reference' / name
            detail = compare(file, reference, case, name)
            detail['sha256'] = hashlib.sha256(file.read_bytes()).hexdigest()
            result['files'][name] = detail
            passed &= detail['passed']
        results['cases'][case] = result
        print(json.dumps({case: result}), flush=True)
    results['passed'] = bool(passed)
    (args.output / 'summary.json').write_text(json.dumps(results, indent=2) + '\n')
    return 0 if passed else 1


if __name__ == '__main__':
    raise SystemExit(main())
