#!/usr/bin/env python3
"""用原 LEO 算例核对生成、加噪复现性，以及 GL/GLML 的实际资料交接。"""
import argparse
import hashlib
import json
import math
from pathlib import Path
import subprocess
import sys
import time


def read_settings(path):
    """读取算例配置的显式关键词，供独立运行目录改写文件路径。"""
    settings = {}

    for line in path.read_text().splitlines():
        if line.strip() and not line.lstrip().startswith('#') and '=' in line:
            key, value = line.split('=', 1)
            settings[key.strip()] = value.strip()

    return settings


def write_settings(path, settings):
    """按原关键词顺序写出本次运行配置。"""
    path.write_text(''.join(f'{key} = {value}\n' for key, value in settings.items()))


def run(command, directory, log):
    """执行真实程序并保存完整输出，失败时给出本次日志位置。"""
    started = time.monotonic()
    result = subprocess.run(command, cwd=directory, capture_output=True, text=True, timeout=120)
    log.write_text(result.stdout + result.stderr)
    assert result.returncode == 0, f'{command[0]} exited {result.returncode}: {log}\n{result.stderr[-2000:]}'
    return time.monotonic() - started


def digest(path):
    """返回文件内容的 SHA-256。"""
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--bin-dir', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True, help='new result directory')
    args = parser.parse_args()
    project = Path(__file__).resolve().parents[1]
    binaries, output = args.bin_dir.resolve(), args.output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    case = project / 'tests/generate_observations'
    summary = {'generator_seconds': {}, 'iod_seconds': {}, 'products': {}, 'noise_products': {}, 'iod_products': {}}
    settings = read_settings(case / 'setting.txt')
    settings.update(file_eph_Keplerian=str(case / 'input/orbit_keplerian.txt'),
                    file_EOP=str(project / 'data/eop/EOP-All_250516.txt'),
                    file_observer=str(case / 'input/observers.txt'),
                    need_output_Cartesian='1', need_output_all='1')

    # 两次使用同一原配置，完整比较 OBS、两类真值和包含不可见点的 all 文件。
    for name in ['theory', 'theory-repeat']:
        folder = output / name
        folder.mkdir()
        config = output / (name + '.txt')
        write_settings(config, dict(settings, dir_output=str(folder)))
        summary['generator_seconds'][name] = run(
            [str(binaries / 'aster_generate_observations'), '-I', str(config)], output,
            output / (name + '-console.txt'))

    theory = output / 'theory'
    products = sorted(theory.iterdir())
    arcs = sorted(theory.glob('obs_*.txt'))
    assert len(arcs) == 6 and sum(len(path.read_text().splitlines()) for path in arcs) == 364
    assert len(next(theory.glob('all_*.txt')).read_text().splitlines()) == 8641

    for path in products:
        assert path.read_bytes() == (output / 'theory-repeat' / path.name).read_bytes(), path
        summary['products'][path.name] = digest(path)

    # 保留原 Python MT19937/gauss 随机流、内容摘要派生种子和两个样本；反转输入顺序复核稳定性。
    for name, sources in [('noise', arcs), ('noise-repeat', list(reversed(arcs)))]:
        run([sys.executable, str(project / 'tools/add_angular_noise.py'), *map(str, sources),
             '--sigma-arcsec', '3', '--samples', '2', '--seed', '42', '--output', str(output / name)],
            output, output / (name + '-console.txt'))

    for path in sorted((output / 'noise').glob('sample_*/*.txt')):
        relative = path.relative_to(output / 'noise')
        assert path.read_bytes() == (output / 'noise-repeat' / relative).read_bytes(), relative
        summary['noise_products'][str(relative)] = digest(path)
        source_rows = [line.split() for line in (theory / path.name).read_text().splitlines()]
        rows = [line.split() for line in path.read_text().splitlines()]
        assert len(rows) == len(source_rows) and all(a[:4] == b[:4] for a, b in zip(rows, source_rows))

    # 复用现有 GL/GLML 的十条观测配置，只改输入、输出路径，并显式给出已有的回退种子。
    summary['iod'] = {}

    for model in ['gl', 'glml']:
        base = read_settings(project / f'tests/iod_{model}/setting.txt')

        for label, source in [('theory', arcs[0]), ('noise', output / 'noise/sample_000001' / arcs[0].name)]:
            name = f'iod-{model}-{label}'
            folder = output / name
            folder.mkdir()
            values = dict(base, file_observation=str(source), file_observer=settings['file_observer'],
                          file_EOP=settings['file_EOP'], random_seed='42', display_iterations='0',
                          file_residual=str(folder / 'residuals.txt'), file_report=str(folder / 'report.txt'))
            values['file_orbit_solution_' + model.upper()] = str(folder / 'solution.txt')

            if model == 'glml':
                values['file_candidates_GLML'] = ''

            config = folder / 'setting.txt'
            write_settings(config, values)
            summary['iod_seconds'][name] = run([str(binaries / f'aster_iod_{model}'), '-I', str(config)],
                                               folder, folder / 'console.txt')
            rows = [line.split() for line in (folder / 'residuals.txt').read_text().splitlines()
                    if line.strip() and not line.lstrip().startswith('#')]
            assert len(rows) == 10 and all(len(row) == 6 for row in rows), name
            original = [line.split() for line in source.read_text().splitlines()[:10]]
            assert all(a[1:4] == b[1:4] and abs(float(a[0])-float(b[0])) < 1e-9
                       and all(math.isfinite(float(value)) for value in a[4:])
                       for a, b in zip(rows, original)), name
            for product in ['solution.txt', 'residuals.txt']:
                summary['iod_products'][name + '/' + product] = digest(folder / product)

            solution = read_settings(folder / 'solution.txt')
            assert solution['object_id'] == 'LEO_2body', solution
            state = [float(value) for value in solution['reference_Cartesian'].split()]
            assert len(state) == 6 and all(math.isfinite(value) for value in state), name
            summary['iod'][name] = {'residual_rows': len(rows),
                                    'solution_flag': solution.get('solution_flag'),
                                    'source_solution_flag': solution.get('source_solution_flag')}

    version = subprocess.run([str(binaries / 'aster_generate_observations'), '--version'],
                             capture_output=True, text=True, check=True).stdout.strip()
    summary.update(version=version, arcs=6, visible_observations=364, all_epoch_rows=8641,
                   noise_sigma_arcsec=3, noise_base_seed=42, noise_samples=2)
    (output / 'summary.json').write_text(json.dumps(summary, ensure_ascii=False, indent=2) + '\n')
    print(json.dumps({key: value for key, value in summary.items()
                      if key not in ['products', 'noise_products', 'iod_products']}, ensure_ascii=False, indent=2))


if __name__ == '__main__':
    main()
