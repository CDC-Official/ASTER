"""检查客户入口的 Cartesian 输出、末步、GM 覆盖和输入错误。"""
import argparse
import csv
import math
from pathlib import Path
import subprocess

parser = argparse.ArgumentParser()
parser.add_argument('--binary', type=Path, required=True)
parser.add_argument('--output', type=Path, required=True)
args = parser.parse_args()
root = Path(__file__).resolve().parents[2]
args.output = args.output.resolve()
args.output.mkdir(parents=True, exist_ok=True)
base = {}
for line in (Path(__file__).parent / 'setting.txt').read_text().splitlines():
    line = line.split('#', 1)[0].strip()
    if line:
        key, value = line.split('=', 1)
        base[key.strip()] = value.strip()
for key in ('file_JPLEPH', 'file_EOP', 'file_lunar_gravity_field'):
    base[key] = str((Path(__file__).parent / base[key]).resolve())


def run(name, changes, failure=None):
    """运行一组显式配置，失败样例同时核对退出码和原因。"""
    directory = args.output / name
    directory.mkdir(exist_ok=True)
    setting = base | {'propagation_length': '125', 'file_eph_Cartesian': str(directory / 'cartesian.csv'),
                      'file_report': str(directory / 'summary.txt'), 'file_programlog': str(directory / 'program.log')} | changes
    path = directory / 'setting.txt'
    path.write_text('\n'.join(f'{key} = {value}' for key, value in setting.items()) + '\n')
    result = subprocess.run([str(args.binary.resolve()), '-I', str(path)], capture_output=True, text=True)
    (directory / 'console.txt').write_text(result.stdout + result.stderr)
    if failure:
        assert result.returncode != 0, name
        assert failure in result.stdout + result.stderr, name
        return None
    assert result.returncode == 0, (name, result.stdout, result.stderr)
    with (directory / 'cartesian.csv').open() as stream:
        rows = [{key: float(value) for key, value in row.items()} for row in csv.DictReader(stream)]
    assert all(math.isfinite(value) for row in rows for value in row.values()), name
    return rows, (directory / 'summary.txt').read_text()


# 非整步时长必须精确到达指定秒数，默认和显式零 EOP 产生同一轨迹。
rows, summary = run('tail_step', {})
assert [row['elapsed_TT_s'] for row in rows] == [0, 60, 120, 125]
assert abs(rows[-1]['MJD_UTC'] - (61892 + 125 / 86400)) < 1e-10
assert 'zero_EOP_queries=39' in summary
zero_rows, _ = run('explicit_zero', {'file_EOP': ''})
assert rows == zero_rows

# 独立圆轨道解析解检查有量纲状态、TT 秒和归一化尺度。
radius = 7000.0
mu = 398600.4415
speed = math.sqrt(mu / radius)
circular, _ = run('circular', {
    'earth_zonal_order': '0', 'earth_tesseral_degree': '0', 'solar_gravity': '0', 'lunar_gravity': '0',
    'reference_Cartesian': f'{radius} 0 0 0 {speed} 0', 'propagation_length': '600', 'integrator_stepsize': '5',
})
phase = math.sqrt(mu / radius**3) * 600
expected = (radius * math.cos(phase), radius * math.sin(phase), 0,
            -speed * math.sin(phase), speed * math.cos(phase), 0)
keys = ('x_km', 'y_km', 'z_km', 'vx_km_s', 'vy_km_s', 'vz_km_s')
assert max(abs(circular[-1][key] - value) for key, value in zip(keys[:3], expected[:3])) < 1e-7
assert max(abs(circular[-1][key] - value) for key, value in zip(keys[3:], expected[3:])) < 1e-10

# 改变月球文件头 GM 而保留 C/S，结果仍由入口显式指定的 JPL GM 决定。
# 这验证实际加速度使用覆盖值，而非只改变摘要打印。
field = Path(base['file_lunar_gravity_field']).read_text()
changed_field = args.output / 'changed_file_gm.txt'
lines = field.splitlines()
header = lines[0].split(',')
assert len(header) > 2, lines[0]
header[1] = '9.80559961386338e12'
lines[0] = ','.join(header)
changed_field.write_text('\n'.join(lines) + '\n')
moon_rows, moon_summary = run('lunar_harmonics', {'use_lunar_NSG': '1', 'lunar_zonal_order': '12', 'lunar_tesseral_degree': '12'})
changed_rows, _ = run('lunar_changed_header', {'use_lunar_NSG': '1', 'lunar_zonal_order': '12', 'lunar_tesseral_degree': '12', 'file_lunar_gravity_field': str(changed_field)})
assert moon_rows == changed_rows
assert 'Moon_field_effective_GM_TT_km3_s2=4902.8001906' in moon_summary
assert moon_rows[-1] != rows[-1]

# 零时长只保存初态；非法设置和错误历表必须明确失败。
initial, _ = run('initial_only', {'propagation_length': '0'})
assert len(initial) == 1
run('bad_step', {'integrator_stepsize': '0'}, 'integrator_stepsize>0')
run('bad_degree', {'earth_zonal_order': '81'}, 'gravity degrees 0..80')
run('bad_ephemeris', {'file_JPLEPH': str(root / 'data/ephemeris/JPLEPH_C.430')}, 'Require native binary DE440')
run('bad_keyword', {'earth_degre': '12'}, 'Unknown keyword')
run('bad_eop_path', {'file_EOP': str(args.output / 'missing_eop.txt')}, 'Failed to open EOP')

# 配置仅使用原主程序已有关键词；初态继续写在 reference_Cartesian 一行中。
original_keys = {line.split('=', 1)[0].strip() for line in (root / 'tests/opn_cartesian/setting.txt').read_text().splitlines() if '=' in line and not line.lstrip().startswith('#')}
assert set(base) <= original_keys | {'file_programlog'}
run('reject_renamed_keyword', {'x_km': '1'}, 'Unknown keyword')
calendar_rows, _ = run('calendar_epoch', {'reference_epoch': '2028 05 01 00 00 00'})
jd_rows, _ = run('jd_epoch', {'reference_epoch': '2461892.5'})
assert calendar_rows == jd_rows == rows
zonal_rows, zonal_summary = run('separate_orders', {'earth_tesseral_degree': '0'})
assert 'Earth_order_degree=12,0' in zonal_summary
assert zonal_rows[-1] != rows[-1]

print('Runtime checks passed: tail step, zero EOP, circular orbit, lunar GM override, original settings and errors.')
