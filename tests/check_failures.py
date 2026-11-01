#!/usr/bin/env python3
"""检查转换输入契约、配置错误和后期写出失败。"""
import argparse
import json
import os
from pathlib import Path
import subprocess
import tempfile
from run_case import prepare


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--bin', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    binaries = args.bin.resolve()
    args.output.mkdir(parents=True, exist_ok=True)
    root = Path(tempfile.mkdtemp(prefix='run-', dir=args.output.resolve()))
    records = []

    def execute(name, target, argv, good=False, cwd=None):
        proc = subprocess.run([str(binaries / target)] + argv, cwd=cwd or root,
                              capture_output=True, env=dict(os.environ, OMP_NUM_THREADS='2'), timeout=30)
        (root / (name + '.console.txt')).write_bytes(proc.stdout + proc.stderr)
        assert (proc.returncode == 0) == good, (name, proc.stderr)
        records.append({'check': name, 'exit': proc.returncode})
        return proc

    tests = Path(__file__).resolve().parent
    for case in ['convert_keplerian_to_cartesian', 'convert_cartesian_to_keplerian']:
        target = 'aster_' + case
        first = next(line for line in (tests / case / 'input.txt').read_text().splitlines() if line.strip())
        unit = ['m'] if case.endswith('cartesian_to_keplerian') else []
        invalid = {'short': ' '.join(first.split()[:6]), 'extra': first + ' 9',
                   'nan': 'nan ' + ' '.join(first.split()[1:]), 'text': first + ' trailing'}
        for label, line in invalid.items():
            path = root / (case + '-' + label + '.txt')
            path.write_text('\n' + line + '\n')
            proc = execute(case + '-' + label, target, [str(path), str(path) + '.out'] + unit)
            assert b'line 2' in proc.stderr
        path = root / (case + '-spaces.txt')
        path.write_text('\n  ' + first + '   \n\n')
        execute(case + '-spaces', target, [str(path), str(path) + '.out'] + unit, good=True)
        assert Path(str(path) + '.out').read_bytes() == (tests / case / 'reference/result.txt').read_bytes()

    cases = json.loads((tests / 'cases.json').read_text())
    for case, info in cases.items():
        proc = execute(case + '-version', info['target'], ['--version'], good=True)
        assert proc.stdout.decode().strip() == info['target'] + ' 1.0.0'
        if not info.get('conversion'):
            execute(case + '-missing-setting', info['target'], ['-I', str(root / 'missing.txt')])

    for label, value in [('negative', '-1'), ('fraction', '1.5'), ('overflow', '2147483648'), ('text', 'abc')]:
        info, _, run, argv = prepare('iod_glml', root / ('seed-' + label))
        with (run / 'setting.txt').open('a') as stream:
            stream.write('\nrandom_seed = ' + value + '\n')
        execute('seed-' + label, info['target'], argv, cwd=run)

    info, _, run, argv = prepare('iod_gl', root / 'late-failure')
    setting = run / 'setting.txt'
    lines = setting.read_text().splitlines()
    setting.write_text('\n'.join('file_report = missing-parent/report.txt' if line.startswith('file_report =') else line for line in lines) + '\n')
    proc = execute('late-report', info['target'], argv, cwd=run)
    assert b'ERROR' in proc.stderr and b'Event Log' not in proc.stderr
    (root / 'checks.json').write_text(json.dumps(records, indent=2))
    print(json.dumps({'checks': len(records), 'directory': str(root)}))


if __name__ == '__main__':
    main()
