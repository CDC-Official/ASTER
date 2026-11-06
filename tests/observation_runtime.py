"""通过真实观测生成入口检查理论观测、真值对齐、输入约束和结果保护。"""
import argparse
import datetime
import math
import re
import struct
from pathlib import Path
import subprocess
import tempfile

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--binary', type=Path, required=True)
parser.add_argument('--project', type=Path, required=True)
parser.add_argument('--jpl', type=Path, help='Optional local JPL binary for numerical-Sun checks')
args = parser.parse_args()
binary = args.binary.resolve()
project = args.project.resolve()
checks = 0


def require(condition, message):
    global checks
    checks += 1
    if not condition:
        raise AssertionError(message)


def values(path):
    result = {}
    for line in path.read_text().splitlines():
        if '=' in line and line.strip() and not line.lstrip().startswith('#'):
            key, value = line.split('=', 1)
            result[key.strip()] = value.strip()
    return result


with tempfile.TemporaryDirectory(prefix='cdc-observation-') as folder:
    root = Path(folder)
    case = project/'tests/generate_observations'
    base = values(case/'setting.txt')
    base.update(file_eph_Keplerian=str(case/'input/orbit_keplerian.txt'),
                file_EOP=str(project/'data/eop/EOP-All_250516.txt'),
                file_observer=str(case/'input/observers.txt'),
                need_output_Keplerian='1', need_output_Cartesian='1')

    def invoke(command, good=True):
        p = subprocess.run([str(binary), *map(str, command)], cwd=root,
                           capture_output=True, text=True, timeout=30)
        log = p.stdout+p.stderr
        require((p.returncode == 0) == good, log[-3000:])
        if not good:
            require('[ERROR]' in log and 'Call Stack' in log, log)
            require('finished successfully' not in log, log)
        return p

    def run(name, changes=None, good=True, out=None):
        d = dict(base)
        for key, value in (changes or {}).items():
            if value is None:
                d.pop(key, None)
            else:
                d[key] = value
        output = out or root/name
        output.mkdir(exist_ok=True)
        d['dir_output'] = str(output)
        config = root/(name+'.txt')
        config.write_text('\n'.join(k+' = '+v for k, v in d.items())+'\n')
        p = invoke(['-I', config], good)
        return output, p

    invoke(['--help'])
    version = invoke(['--version'])
    expected_version = re.search(r'project\(ASTER VERSION ([0-9.]+)', (project/'CMakeLists.txt').read_text()).group(1)
    expected_version += re.search(r'set\(ASTER_VERSION_SUFFIX "([^"]*)"\)',
                                  (project/'CMakeLists.txt').read_text()).group(1)
    require(version.stdout.startswith('aster_generate_observations '+expected_version+'\n')
            and 'backend=' in version.stdout, 'Missing program version or math backend')
    no_arguments = subprocess.run([str(binary)], cwd=root, capture_output=True, text=True, timeout=10)
    require(no_arguments.returncode == 2 and 'Usage:' in no_arguments.stderr,
            'Missing explicit configuration must show usage')

    for command in [['--unknown'], ['-I'], ['-I', root/'missing'], ['--output-dir', root/'invalid']]:
        invoke(command, False)

    # 输出目录由调用方准备；缺少目录时不能偷偷创建。
    missing_dir_setting = root/'missing-output-dir.txt'
    missing_dir_values = dict(base, dir_output=str(root/'not-prepared'))
    missing_dir_setting.write_text('\n'.join(k+' = '+v for k,v in missing_dir_values.items())+'\n')
    invoke(['-I', missing_dir_setting], False)
    require(not (root/'not-prepared').exists(), 'Program created an output directory')

    first, p = run('first')
    second, _ = run('second')
    for flag in ['-input', '--input']:
        invoke([flag, root/'second.txt'])
    files = sorted(f.name for f in first.glob('obs_*.txt'))
    require(bool(files), 'Fixture produced no observations')
    products = [p for p in first.iterdir() if p.name.startswith(('obs_', 'kep_', 'cart_'))]
    require(all(p.read_bytes() == (second/p.name).read_bytes() for p in products), 'Theoretical output changed')
    lines = p.stdout.splitlines()
    require(int(lines[0]) == len(files) and len(lines) == len(files)+1, 'Stdout file list changed')
    require(all(Path(name).is_file() for name in lines[1:]), 'File-list entry is missing')
    for name in files:
        obs = [row.split() for row in (first/name).read_text().splitlines()]
        kep = [[float(v) for v in row.split()] for row in (first/name.replace('obs_', 'kep_', 1)).read_text().splitlines()]
        cart = [[float(v) for v in row.split()] for row in (first/name.replace('obs_', 'cart_', 1)).read_text().splitlines()]
        require(len(obs) == len(kep) == len(cart), 'Truth/observation row counts differ')
        require(all(len(o)==6 and len(k)==7 and len(c)==7
                    and float(format(float(o[0]), '.15g'))==k[0]==c[0]
                    and o[1]=='1' and o[2]==base['object_id'] and o[3]=='1'
                    and 0<=float(o[4])<360 and -90<=float(o[5])<=90
                    and all(math.isfinite(float(v)) for v in [o[0],o[4],o[5]])
                    and all(math.isfinite(v) for v in k+c)
                    and 6000<math.sqrt(sum(v*v for v in c[1:4]))<8000
                    and 5<math.sqrt(sum(v*v for v in c[4:7]))<10
                    for o, k, c in zip(obs, kep, cart)), 'Wrong columns, units or epoch alignment')
    itrf = root/'itrf-stations.txt'
    # WGS84：纬度 41.725892°、经度 86.174633°、高 943 m 对应的独立 ECEF 坐标。
    itrf.write_text('1 2 318101.62456398865 4757397.466240386 4223557.521084583\n')
    alternate, _ = run('itrf-site', {'file_observer':str(itrf)})
    require(all(p.read_bytes()==(alternate/p.name).read_bytes() for p in products), 'ITRS/geodetic station disagreement')
    named = root/'named-stations.txt'
    named.write_text('OPTICAL_A 1 41.725892 86.174633 943.0\n')
    named_run, _ = run('named-site', {'file_observer':str(named), 'observer_id':'OPTICAL_A'})
    require(all(line.split()[3]=='OPTICAL_A' for f in named_run.glob('obs_*.txt') for line in f.read_text().splitlines()),
            'String observer ID was changed or padded in OBS')
    # 完整文件一一对应输入，前六列须复用弧段结果；开关不能改变任何旧产品。
    input_records = [row.split() for row in (case/'input/orbit_keplerian.txt').read_text().splitlines() if row.strip()]

    def complete(folder, expected_input=input_records):
        paths = list(folder.glob('all_*.txt'))
        require(len(paths)==1, 'Expected one all_ file')
        text = paths[0].read_text().splitlines()
        fields = [line.split() for line in text]
        require(len(fields)==len(expected_input) and all(len(row)==8 for row in fields), 'Complete-file row/column count changed')
        require(all(float(row[0])==float(source[1]) and int(row[7])==int(source[0])
                    for row, source in zip(fields, expected_input)), 'Input epoch/ID order changed')
        require(all(row[1]=='1' and row[2]==base['object_id'] and row[3]=='1' and row[6] in ['0','1']
                    and all(math.isfinite(float(v)) for v in [row[0],*row[4:]])
                    and 0<=float(row[4])<360 and -90<=float(row[5])<=90
                    for row in fields),
                'Malformed complete-file values')
        arc_rows = {float(line.split()[0]):line for file in folder.glob('obs_*.txt') for line in file.read_text().splitlines()}
        visible = {float(row[0]) for row in fields if row[6]=='1'}
        require(visible==arc_rows.keys(), 'Complete visibility flags disagree with arcs')
        require(all(line.startswith(arc_rows[float(row[0])]) for line, row in zip(text, fields) if row[6]=='1'),
                'Visible row prefix differs from the arc file')
        return paths[0], fields

    require(not list(first.glob('all_*.txt')), 'Default unexpectedly creates a complete file')
    missing_switch, _ = run('missing-all-switch', {'need_output_all':None})
    require(not list(missing_switch.glob('all_*.txt')) and
            all(p.read_bytes()==(missing_switch/p.name).read_bytes() for p in products), 'Missing switch broke compatibility')
    full, full_process = run('complete', {'need_output_all':'1'})
    full_file, full_rows = complete(full)
    require(all(p.read_bytes()==(full/p.name).read_bytes() for p in products), 'Complete output changed old arc/truth products')
    require(len(full_process.stdout.splitlines())==len(files)+1 and 'all_' not in full_process.stdout,
            'Complete file leaked into the stdout arc list')
    require(str(full_file) in full_process.stderr, 'Complete path is absent from the run summary')
    repeated_full, _ = run('complete-repeat', {'need_output_all':'1'})
    repeated_file, _ = complete(repeated_full)
    require(full_file.read_bytes()==repeated_file.read_bytes(), 'Complete output is not reproducible')
    empty_sky, empty_process = run('complete-none-visible', {'need_output_all':'1', 'use_station_night':'1', 'max_sun_elevation_deg':'-90'})
    _, empty_rows = complete(empty_sky)
    require(all(row[6]=='0' for row in empty_rows) and empty_process.stdout=='0\n', 'All-invisible output is incomplete')
    filtered, _ = run('filtered-without-complete', {'use_station_night':'1', 'use_earth_shadow':'1'})
    filtered_full, _ = run('filtered-with-complete', {'use_station_night':'1', 'use_earth_shadow':'1', 'need_output_all':'1'})
    complete(filtered_full)
    require(all(p.read_bytes()==(filtered_full/p.name).read_bytes() for p in filtered.iterdir()
                if p.name.startswith(('obs_', 'kep_', 'cart_'))), 'Complete output changes filtered theoretical arcs')
    sparse_input = root/'sparse-id-input.txt'
    sparse_rows = [list(row) for row in input_records[:3]]
    for row, value in zip(sparse_rows, ['17','25','100']): row[0]=value
    sparse_input.write_text('\n\n'.join(' '.join(row) for row in sparse_rows)+'\n')
    sparse, _ = run('complete-sparse-id', {'file_eph_Keplerian':str(sparse_input), 'need_output_all':'1'})
    complete(sparse, sparse_rows)

    sentinel = first/'keep.txt'
    sentinel.write_text('preserve')
    run('repeat', out=first)
    require(sentinel.read_text()=='preserve' and not (root/'first-rerun-1').exists(),
            'Program changed the explicit output directory')
    require(not (first/'config.original.txt').exists() and not (first/'paths.txt').exists(),
            'Unexpected automatic run snapshots')

    for i, changes in enumerate([
        {'file_eph_Cartesian':base['file_eph_Keplerian']},
        {'file_eph_Keplerian':''}, {'observation_type':'2'},
        {'min_elevation_deg':'90'},
        {'random_seed':'123'}, {'noise_sigma_as':'1'},
        {'observation_type':None},
        {'need_output_Cartesian':'2'}, {'need_output_all':'2'}, {'observer_id':'missing'},
        {'is_input_normalized':'2'}, {'object_id':'../target'},
    ]):
        run('bad-'+str(i), changes, False)
    # 恢复原读取语义后，缺项或不能解析的高度角警告并读为零；零是有效高度角。
    for name, value in [('missing-elevation', None), ('invalid-elevation', 'nan')]:
        output, result = run(name, {'min_elevation_deg':value})
        require('[WARNING]' in result.stdout + result.stderr, 'Missing or invalid elevation must warn')
    empty = root/'empty-ephemeris.txt'
    empty.write_text('')
    run('empty', {'file_eph_Keplerian':str(empty)}, False)
    source = (case/'input/orbit_keplerian.txt').read_text().splitlines()
    malformed = root/'malformed-ephemeris.txt'
    malformed.write_text(source[0]+' extra\n')
    run('malformed', {'file_eph_Keplerian':str(malformed)}, False)
    malformed.write_text(source[1]+'\n'+source[0]+'\n')
    run('decreasing', {'file_eph_Keplerian':str(malformed)}, False)
    # 首历元有覆盖、后续没有覆盖时也必须失败，防止恢复为整段复用首值。
    day = math.floor(float(source[0].split()[1]))
    limited = root/'limited-eop-data.txt'
    rows = []
    for mjd in [day, day+1]:
        date = datetime.date(1858, 11, 17)+datetime.timedelta(days=mjd)
        rows.append(f'{date.year} {date.month} {date.day} {mjd} 0.1 0.2 0.3 0 0 0 0 0 37')
    limited.write_text('BEGIN OBSERVED\n'+'\n'.join(rows)+'\nEND OBSERVED\n')
    _, limited_result = run('limited-eop', {'file_EOP':str(limited)}, False)
    require('Invalid EOP coverage' in limited_result.stderr, 'Wrong EOP failure: '+limited_result.stderr)

    # 关闭筛选时无需太阳文件；打开后检验两个条件的交集与门限单调性。
    disabled, _ = run('unused-jpl', {'solar_ephemeris':'1', 'file_JPLEPH':'missing-jpl'})
    require(all((first/f).read_bytes()==(disabled/f).read_bytes() for f in files), 'Unused Sun settings changed observations')

    def epochs(folder):
        return {float(row.split()[0]) for f in folder.glob('obs_*.txt') for row in f.read_text().splitlines()}

    def check_filters(mode, jpl=None):
        sets = {}
        for night, shadow in [(0, 0), (1, 0), (0, 1), (1, 1)]:
            changes = {'solar_ephemeris':str(mode), 'use_station_night':str(night),
                       'use_earth_shadow':str(shadow), 'max_sun_elevation_deg':'-6', 'need_output_all':'1'}
            if jpl: changes['file_JPLEPH'] = str(jpl)
            output, _ = run(f'sun-{mode}-{night}-{shadow}', changes)
            complete(output)
            sets[night, shadow] = epochs(output)
        require(sets[1,0] <= sets[0,0] and sets[0,1] <= sets[0,0], 'Filter added observations')
        require(sets[1,1] == sets[1,0] & sets[0,1], 'Combined filters differ from the intersection')
        return sets

    analytic = check_filters(2)
    darker, _ = run('darker', {'use_station_night':'1', 'max_sun_elevation_deg':'-18'})
    require(epochs(darker) <= analytic[1,0], 'Stricter darkness threshold added observations')
    higher, _ = run('higher', {'min_elevation_deg':'10'})
    require(epochs(higher) <= analytic[0,0], 'Higher elevation threshold added observations')
    for i, changes in enumerate([
        {'use_station_night':'2'}, {'use_earth_shadow':'-1'}, {'use_earth_shadow':'1', 'solar_ephemeris':'0'},
        {'use_station_night':'1', 'max_sun_elevation_deg':None},
        {'use_station_night':'1', 'max_sun_elevation_deg':'0'},
        {'use_earth_shadow':'1', 'solar_ephemeris':'1', 'file_JPLEPH':''},
        {'use_earth_shadow':'1', 'solar_ephemeris':'1', 'file_JPLEPH':'missing'},
    ]):
        run('bad-filter-'+str(i), changes, False)
    if args.jpl:
        numerical = check_filters(1, args.jpl.resolve())
        # 保留真实文件的合法头部，构造观测历元在覆盖外的短文件。
        # 仅检验 app 在调用底层插值器前检查范围，不用于轨道精度验证。
        clipped = root/'out-of-range-jpl-data.bin'
        with args.jpl.open('rb') as stream:
            header = bytearray(stream.read(3*8144))
        offset = 3*84+400*6  # 现有 packed rec1 中 ss[3] 的字节偏移。
        day = float(source[0].split()[1])+2400000.5
        struct.pack_into('=ddd', header, offset, day+100, day+132, 32.0)
        clipped.write_bytes(header)
        _, invalid = run('out-of-range-jpl', {'solar_ephemeris':'1', 'use_earth_shadow':'1', 'file_JPLEPH':str(clipped)}, False)
        require('outside JPL coverage' in invalid.stderr and not invalid.stdout, 'JPL coverage failure did not reach tracer cleanly')
        print('Optical filter counts [off, night, shadow, both]:',
              {label:[len(sets[k]) for k in [(0,0),(1,0),(0,1),(1,1)]]
               for label, sets in [('analytical', analytic), ('JPL', numerical)]})

    # 从已写出的真值构造同一段的两种输入，核对真实入口的单位转换。
    name = files[0]
    kep_rows = (first/name.replace('obs_', 'kep_', 1)).read_text().splitlines()
    cart_rows = (first/name.replace('obs_', 'cart_', 1)).read_text().splitlines()
    for kind, rows in [('kep', kep_rows), ('cart', cart_rows)]:
        (root/(kind+'.txt')).write_text('\n'.join(str(i+1)+' '+row for i, row in enumerate(rows))+'\n')
    kep_run, _ = run('kep-roundtrip', {'file_eph_Keplerian':str(root/'kep.txt')})
    cart_run, _ = run('cart-roundtrip', {'file_eph_Keplerian':'', 'file_eph_Cartesian':str(root/'cart.txt')})
    a = [row.split() for row in next(kep_run.glob('obs_*.txt')).read_text().splitlines()]
    b = [row.split() for row in next(cart_run.glob('obs_*.txt')).read_text().splitlines()]
    require(len(a)==len(b) and all(abs(float(x[4])-float(y[4]))<5e-8 and abs(float(x[5])-float(y[5]))<5e-8 for x,y in zip(a,b)), 'Cartesian/Keplerian input disagreement')

print(f'Observation runtime checks passed: {checks}')
