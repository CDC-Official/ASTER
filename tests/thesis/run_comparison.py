#!/usr/bin/env python3
"""对已提交的 v1.4.0 与论文实验分支做串行交错计时和逐字节科学结果比较。"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import statistics
import subprocess
import time


def sha256(path):
    """返回给定文件的 SHA-256，用于固定科学输入和输出。"""
    return hashlib.sha256(path.read_bytes()).hexdigest()


def git(source, *args):
    """读取指定源码目录的 Git 信息。"""
    return subprocess.check_output(['git', '-C', str(source), *args], text=True).strip()


def prepare_setting(fixture, run, overrides):
    """解析已有 setting，固定输入绝对路径并保存本次实际配置。"""
    lines = []
    inputs = {}
    for line in (fixture / 'setting.txt').read_text().splitlines():
        if '=' not in line or line.lstrip().startswith('#'):
            continue
        key, value = map(str.strip, line.split('=', 1))
        if key in ('file_observation', 'file_observer', 'file_EOP'):
            path = (fixture / value).resolve()
            value = str(path)
            inputs[key] = {'path': value, 'sha256': sha256(path)}
        value = overrides.get(key, value)
        lines.append(f'{key} = {value}')
    if 'random_seed' not in [line.split('=')[0].strip() for line in lines]:
        lines.append('random_seed = 0')
    setting = run / 'setting.txt'
    setting.write_text('\n'.join(lines) + '\n')
    return setting, inputs


def run_once(binary, fixture, run, overrides, threads):
    """运行独立进程，保存配置、控制台、精确墙钟、资源和输出校验。"""
    run.mkdir(parents=True, exist_ok=False)
    setting, inputs = prepare_setting(fixture, run, overrides)
    env = dict(os.environ, OMP_NUM_THREADS=str(threads), OMP_DYNAMIC='FALSE',
               OMP_PROC_BIND='spread', OMP_PLACES='cores')
    command = [str(binary), '-I', str(setting)]
    if Path('/usr/bin/time').exists() and os.uname().sysname == 'Linux':
        command = ['/usr/bin/time', '-f', '%e %U %S %M', '-o', str(run / 'process_time.txt'), *command]
    started = time.perf_counter_ns()
    with (run / 'console.txt').open('wb') as console:
        result = subprocess.run(command, cwd=run, env=env, stdout=console, stderr=subprocess.STDOUT)
    wall_s = (time.perf_counter_ns() - started) / 1e9
    if result.returncode:
        raise RuntimeError(f'GLML failed in {run}, exit={result.returncode}')
    names = ('glml_candidates.csv', 'glml_solution.txt', 'glml_residuals.txt')
    hashes = {name: sha256(run / name) for name in names}
    report = (run / 'glml_report.txt').read_text()
    # 路径、事件行号、耗时和内存随运行变化；迭代文字单独比较。
    iterations = re.findall(r'\n-+ Iteration \d+ -+\n.*?(?=\n-+ Iteration \d+ -+\n|\n=+|\Z)', report, re.S)
    iteration_text = ''.join(iterations)
    (run / 'iteration_report.txt').write_text(iteration_text)
    selected = re.search(r'^selected_candidate\s*=\s*(\d+)', (run / 'glml_solution.txt').read_text(), re.M)
    data = {'directory': str(run), 'wall_s': wall_s, 'hashes': hashes, 'inputs': inputs,
            'iteration_report_sha256': sha256(run / 'iteration_report.txt'),
            'report_iteration_count': len(iterations), 'selected_candidate': int(selected[1]) if selected else None}
    console = (run / 'console.txt').read_text()
    generation = re.findall(r'Progress:.*?elapsed ([\d.]+) s', console)
    data['generation_display_s'] = float(generation[-1]) if generation else None
    if (run / 'process_time.txt').exists():
        values = (run / 'process_time.txt').read_text().split()
        data['gnu_time'] = dict(zip(('wall_s', 'user_s', 'system_s', 'max_rss_kib'), map(float, values)))
    (run / 'run.json').write_text(json.dumps(data, indent=2) + '\n')
    return data


def main():
    """固定双方版本后依次运行三组输入；预热单独保存，计时不与其他组并发。"""
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('baseline-source', 'candidate-source', 'baseline-bin', 'candidate-bin', 'output'):
        parser.add_argument('--' + name, type=Path, required=True)
    parser.add_argument('--threads', type=int, nargs='+', default=[32, 64])
    parser.add_argument('--repeats', type=int, default=5)
    parser.add_argument('--cases', nargs='+', choices=['sensitive10', 'short20', 'long100'],
                        default=['sensitive10', 'short20', 'long100'])
    args = parser.parse_args()
    if args.repeats < 1 or any(n < 1 for n in args.threads):
        parser.error('repeats and threads must be positive')
    args.output = args.output.resolve()
    args.output.mkdir(parents=True, exist_ok=False)
    sources = {'baseline': args.baseline_source.resolve(), 'candidate': args.candidate_source.resolve()}
    binaries = {'baseline': (args.baseline_bin / 'aster_iod_glml').resolve(),
                'candidate': (args.candidate_bin / 'aster_iod_glml').resolve()}
    metadata = {}
    for mode, source in sources.items():
        if git(source, 'diff', '--name-only', 'HEAD'):
            raise RuntimeError(f'{mode} has uncommitted tracked changes')
        metadata[mode] = {'source': str(source), 'commit': git(source, 'rev-parse', 'HEAD'),
                          'binary': str(binaries[mode]), 'binary_sha256': sha256(binaries[mode]),
                          'version': subprocess.check_output([str(binaries[mode]), '--version'], text=True).strip()}
    if metadata['baseline']['commit'] != git(sources['baseline'], 'rev-parse', 'v1.4.0^{}'):
        raise RuntimeError('Baseline must be the published v1.4.0 commit')
    fixture_root = sources['baseline'] / 'tests'
    fixtures = {'sensitive10': (fixture_root / 'iod_glml', {'display_iterations': '0'}),
                'short20': (fixture_root / 'math_backend_cases/cases/e070_s3_n020_m0002', {}),
                'long100': (fixture_root / 'math_backend_cases/cases/e010_s0.5_n100_m0001', {})}
    summary = {'metadata': metadata, 'uname': list(os.uname()), 'threads': args.threads,
               'repeats': args.repeats, 'timing': 'independent processes; one excluded warm-up per version; alternating baseline/candidate order; perf_counter_ns includes process launch and exit',
               'groups': [], 'all_equal': True}
    for case in args.cases:
        fixture, overrides = fixtures[case]
        for threads in args.threads:
            group = {'case': case, 'threads': threads, 'warmup': {}, 'runs': {'baseline': [], 'candidate': []}}
            for mode in sources:
                group['warmup'][mode] = run_once(binaries[mode], fixture, args.output / case / f't{threads}' / f'{mode}-warmup', overrides, threads)
            for repeat in range(args.repeats):
                modes = ['baseline', 'candidate'] if repeat % 2 == 0 else ['candidate', 'baseline']
                for mode in modes:
                    row = run_once(binaries[mode], fixture, args.output / case / f't{threads}' / f'{mode}-{repeat:02}', overrides, threads)
                    group['runs'][mode].append(row)
                    print(json.dumps({'case': case, 'threads': threads, 'repeat': repeat, 'mode': mode,
                                      'wall_s': row['wall_s'], 'selected': row['selected_candidate']}), flush=True)
            first = group['warmup']['baseline']
            group['scientific_files_equal'] = all(row['hashes'] == first['hashes']
                for mode in sources for row in [group['warmup'][mode], *group['runs'][mode]])
            group['iteration_reports_equal'] = all(row['iteration_report_sha256'] == first['iteration_report_sha256'] and row['report_iteration_count'] > 0
                for mode in sources for row in [group['warmup'][mode], *group['runs'][mode]])
            group['median_s'] = {mode: statistics.median(row['wall_s'] for row in group['runs'][mode]) for mode in sources}
            group['speedup'] = group['median_s']['baseline'] / group['median_s']['candidate']
            summary['all_equal'] &= group['scientific_files_equal'] and group['iteration_reports_equal']
            summary['groups'].append(group)
            (args.output / 'summary.json').write_text(json.dumps(summary, indent=2) + '\n')
            print('GROUP ' + json.dumps({k: v for k, v in group.items() if k not in ('runs', 'warmup')}), flush=True)
    if not summary['all_equal']:
        raise SystemExit('Scientific outputs or numeric iteration report differed; inspect summary.json')


if __name__ == '__main__':
    main()
