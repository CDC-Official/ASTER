#!/usr/bin/env python3
"""运行全部主程序的代表算例，保存科学输出校验值，供 Mac/Linux 发布对照。"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import platform
import re
import subprocess
from run_case import prepare

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--bin-dir', type=Path, required=True)
parser.add_argument('--output', type=Path, required=True)
args = parser.parse_args()
root = Path(__file__).resolve().parent.parent
binary = args.bin_dir.resolve()
output = args.output.resolve()
output.mkdir(parents=True, exist_ok=False)
version = re.search(r'project\(ASTER VERSION ([0-9.]+)', (root / 'CMakeLists.txt').read_text()).group(1)
manifest = {'platform': platform.platform(), 'version': version,
            'commit': subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=root, text=True).strip(),
            'threads': 2, 'programs': {}, 'scientific_files': {}}

def run_program(name, command, directory):
    """运行一次完整主程序，保留日志并检查退出状态。"""
    run = subprocess.run([str(binary / name), *map(str, command)], cwd=directory,
                         env=dict(os.environ, OMP_NUM_THREADS='2', OMP_DYNAMIC='FALSE'),
                         capture_output=True, timeout=180)
    (directory / 'console.txt').write_bytes(run.stdout + run.stderr)
    if run.returncode != 0:
        raise RuntimeError(f'{name} failed; see {directory / "console.txt"}')

def record_file(case, path):
    """记录供两端逐字节比较的科学文件，日志和含机器路径的报告另外保存。"""
    data = path.read_bytes()
    manifest['scientific_files'][case + '/' + path.name] = {
        'sha256': hashlib.sha256(data).hexdigest(), 'bytes': len(data),
        'rows': len(data.splitlines()), 'path': str(path)}

cases = json.loads((root / 'tests/cases.json').read_text())
programs = [info['target'] for info in cases.values()] + ['aster_generate_observations', 'aster_generate_cislunar_observations', 'aster_convert_coordinates', 'aster_convert_states']

# 每个入口检查版本和帮助，再运行已有的小算例。
for name in programs:
    v = subprocess.run([str(binary / name), '--version'], capture_output=True, text=True, check=True)
    h = subprocess.run([str(binary / name), '--help'], capture_output=True, text=True, check=True)
    assert v.stdout.splitlines()[0] == f'{name} {version}', (name, v.stdout)
    assert name in h.stdout + h.stderr, name
    manifest['programs'][name] = {'version_output': v.stdout.strip(), 'help_passed': True}

for case, _ in cases.items():
    info, fixture, directory, command = prepare(case, output / case)
    run_program(info['target'], command, directory)

    for name in info['references']:
        record_file(case, directory / name)

# 地基测角同时保存各弧的观测与两种轨道真值，输出由同一配置准备。
case = root / 'tests/generate_observations'
directory = output / 'generate_observations'
directory.mkdir()
settings = {k.strip(): v.strip() for line in (case / 'setting.txt').read_text().splitlines()
            if '=' in line and not line.lstrip().startswith('#') for k, v in [line.split('=', 1)]}
for key in ['file_eph_Cartesian', 'file_eph_Keplerian', 'file_EOP', 'file_observer', 'file_JPLEPH']:
    if settings.get(key): settings[key] = str((case / settings[key]).resolve())
settings.update(dir_output=str(directory), need_output_Keplerian='1', need_output_Cartesian='1', need_output_all='1')
config = directory / 'setting.txt'
config.write_text('\n'.join(f'{key} = {value}' for key, value in settings.items()) + '\n')
run_program('aster_generate_observations', ['-I', config], directory)
for path in sorted(directory.iterdir()):
    if path.name.startswith(('obs_', 'kep_', 'cart_', 'all_')): record_file('generate_observations', path)

# 两个方向各生成一条测距；采用零噪声，比较观测及全部物理分项。
case = root / 'tests/generate_cislunar_observations'
directory = output / 'generate_cislunar_observations'
directory.mkdir()
settings = {k.strip(): v.strip() for line in (case / 'aster.settings').read_text().splitlines()
            if '=' in line for k, v in [line.split('=', 1)]}
for key in ['file_navigation_state', 'file_user_state', 'file_receive_events', 'file_JPLEPH']:
    settings[key] = str((case / settings[key]).resolve())
settings.update(directory_output=str(directory / 'output'), noise_sigma_m='0')
config = directory / 'setting.txt'
config.write_text('\n'.join(f'{key} = {value}' for key, value in settings.items()) + '\n')
run_program('aster_generate_cislunar_observations', ['-I', config], directory)
for name in ['observations.txt', 'components.csv']: record_file('generate_cislunar_observations', directory / 'output' / name)

# 六个测站转换方向复用同一小算例，保留全部四列科学结果。
directory = output / 'convert_coordinates'
directory.mkdir()
source_paths = {'geodetic': root / 'tests/convert_coordinates/geodetic.txt',
                'itrs': directory / 'geodetic_to_itrs.txt',
                'gcrs': directory / 'geodetic_to_gcrs.txt'}
for source, target in [('geodetic', 'itrs'), ('geodetic', 'gcrs'),
                       ('itrs', 'geodetic'), ('itrs', 'gcrs'),
                       ('gcrs', 'itrs'), ('gcrs', 'geodetic')]:
    path = directory / (source + '_to_' + target + '.txt')
    command = [source, target, source_paths[source], path]
    if 'gcrs' in (source, target): command.append(root / 'data/eop/EOP-All_250516.txt')
    run_program('aster_convert_coordinates', command, directory)
    record_file('convert_coordinates', path)

# 三种轨道状态的六个方向使用两条现有小输入。
directory = output / 'convert_states'
directory.mkdir()
source_paths = {'kep': root / 'tests/convert_states/keplerian_km.txt',
                'cart': directory / 'kep_to_cart.txt',
                'nsg2': directory / 'kep_to_nsg2.txt'}
for source, target in [('kep', 'cart'), ('kep', 'nsg2'),
                       ('cart', 'kep'), ('cart', 'nsg2'),
                       ('nsg2', 'cart'), ('nsg2', 'kep')]:
    path = directory / (source + '_to_' + target + '.txt')
    run_program('aster_convert_states', [source, target, source_paths[source], path, '--units', 'deg-km'], directory)
    record_file('convert_states', path)

manifest['data_manifest_sha256'] = hashlib.sha256((root / 'data/manifest.json').read_bytes()).hexdigest()
(output / 'summary.json').write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + '\n')
print(f'Passed: {len(programs)} programs, {len(manifest["scientific_files"])} scientific files; {output / "summary.json"}')
