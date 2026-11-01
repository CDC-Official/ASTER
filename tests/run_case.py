#!/usr/bin/env python3
"""在独立目录运行最小算例，并与迁移前保存的科学结果比较。"""
import argparse
import json
import os
from pathlib import Path
import subprocess
import tempfile


def prepare(case, output):
    tests = Path(__file__).resolve().parent
    info = json.loads((tests / 'cases.json').read_text())[case]
    fixture = tests / case
    output.mkdir(parents=True, exist_ok=True)
    run = Path(tempfile.mkdtemp(prefix='run-', dir=output))
    if info.get('conversion'):
        args = [str(fixture / 'input.txt'), str(run / 'result.txt')]
        if info['unit']:
            args.append(info['unit'])
    else:
        lines = []
        for line in (fixture / 'setting.txt').read_text().splitlines():
            if '=' not in line or line.lstrip().startswith('#'):
                continue
            key, value = map(str.strip, line.split('=', 1))
            if value and key in info['input_keys']:
                value = str((fixture / value).resolve())
            lines.append(key + ' = ' + value)
        (run / 'setting.txt').write_text('\n'.join(lines) + '\n')
        args = ['-I', str(run / 'setting.txt')]
    return info, fixture, run, args


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--case', required=True)
    parser.add_argument('--binary', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    info, fixture, run, arguments = prepare(args.case, args.output.resolve())
    process = subprocess.run([str(args.binary.resolve())] + arguments, cwd=run,
                             env=dict(os.environ, OMP_NUM_THREADS='2'), capture_output=True, timeout=110)
    (run / 'console.txt').write_bytes(process.stdout + process.stderr)
    assert process.returncode == 0, (run, process.returncode, process.stderr[-2500:])
    for name in info['references']:
        assert (run / name).read_bytes() == (fixture / 'reference' / name).read_bytes(), (run, name)
    print(json.dumps({'case': args.case, 'directory': str(run), 'files_equal': info['references']}))


if __name__ == '__main__':
    main()
