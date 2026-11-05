#!/usr/bin/env python3
"""检查转换输入契约、配置错误和后期写出失败。"""
import argparse
import json
import os
import re
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
        plain = root / (case + '-plain.txt')
        plain.write_text(first + '\n')
        execute(case + '-plain', target, [str(plain), str(plain) + '.out'] + unit, good=True)
        assert Path(str(path) + '.out').read_bytes() == Path(str(plain) + '.out').read_bytes()

    version_match = re.search(r'project\(ASTER VERSION (\d+\.\d+\.\d+) LANGUAGES',
                              (tests.parent / 'CMakeLists.txt').read_text())
    assert version_match, 'Missing ASTER project version'
    version = version_match.group(1)
    cases = json.loads((tests / 'cases.json').read_text())
    for case, info in cases.items():
        proc = execute(case + '-version', info['target'], ['--version'], good=True)
        version_lines = proc.stdout.decode().splitlines()
        assert version_lines[0] == info['target'] + ' ' + version
        assert len(version_lines) >= 2 and version_lines[1].startswith('backend=')
        if not info.get('conversion'):
            execute(case + '-missing-setting', info['target'], ['-I', str(root / 'missing.txt')])

    for label, value in [('negative', '-1'), ('fraction', '1.5'), ('overflow', '2147483648'), ('text', 'abc')]:
        info, _, run, argv = prepare('iod_glml', root / ('seed-' + label))
        with (run / 'setting.txt').open('a') as stream:
            stream.write('\nrandom_seed = ' + value + '\n')
        execute('seed-' + label, info['target'], argv, cwd=run)

    # 两套地球入口都应通过 Tracer 拒绝错误开关和缺失的配套模型资料。
    tide_failures = {
        'bad-value': {'use_solidtide_perturb': '2'},
        'fraction': {'use_solidtide_perturb': '0.5'},
        'missing-jpl': {'use_solidtide_perturb': '1', 'file_JPLEPH': ''},
        'unreadable-jpl': {'use_solidtide_perturb': '1', 'file_JPLEPH': 'missing-jpl.bin'},
        'missing-c20': {'use_solidtide_perturb': '1', 'earth_zonal_order': '0'},
    }
    for case in ['opn_cartesian', 'opn_nonsingular2']:
        for label, changes in tide_failures.items():
            info, _, run, argv = prepare(case, root / (case + '-tide-' + label))
            setting = run / 'setting.txt'
            values = dict(line.split(' = ', 1) for line in setting.read_text().splitlines())
            values.update(changes)
            setting.write_text(''.join(key + ' = ' + value + '\n' for key, value in values.items()))
            proc = execute(case + '-tide-' + label, info['target'], argv, cwd=run)
            assert b'ERROR' in proc.stderr and b'main' in proc.stderr, (case, label, proc.stderr)

    # 仅开启固体潮时，入口仍须加载 JPL；静态场取 2/0 即可提供所需 C20。
    for case in ['opn_cartesian', 'opn_nonsingular2']:
        info, _, run, argv = prepare(case, root / (case + '-tide-only'))
        setting = run / 'setting.txt'
        values = dict(line.split(' = ', 1) for line in setting.read_text().splitlines())
        values.update(use_solidtide_perturb='1', propagation_length='10', earth_zonal_order='2',
                      earth_tesseral_degree='0', solar_gravity='0', lunar_gravity='0',
                      use_lunar_NSG='0', use_drag_perturb='0', use_SRP_perturb='0')
        setting.write_text(''.join(key + ' = ' + value + '\n' for key, value in values.items()))
        proc = execute(case + '-tide-only', info['target'], argv, good=True, cwd=run)
        assert b'Degree-2 elastic, Sun and Moon' in proc.stderr

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
