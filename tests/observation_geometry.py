"""用教材原始角度公式及门限两侧的人工几何，核对真实 app 的筛选结果。"""
import argparse
import json
import math
from pathlib import Path
import subprocess
import tempfile

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--binary', type=Path, required=True)
parser.add_argument('--probe', type=Path, required=True)
parser.add_argument('--project', type=Path, required=True)
parser.add_argument('--jpl', type=Path)
args = parser.parse_args()
project, binary, probe = args.project.resolve(), args.binary.resolve(), args.probe.resolve()
case = project/'tests/generate_observations'
eop = project/'data/eop/EOP-All_250516.txt'
checks = 0


def require(condition, message):
    global checks
    checks += 1
    if not condition:
        raise AssertionError(message)


def norm(v):
    return math.sqrt(sum(x*x for x in v))


def unit(v):
    n = norm(v)
    return [x/n for x in v]


def dot(a, b):
    return sum(x*y for x, y in zip(a, b))


def cross(a, b):
    return [a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]]


def angle(a, b):
    return math.acos(max(-1.0, min(1.0, dot(a, b)/norm(a)/norm(b))))


with tempfile.TemporaryDirectory(prefix='cdc-optical-geometry-') as folder:
    root = Path(folder)
    summaries = {}
    all_probe_rows = {}
    for mode in [2]+([1] if args.jpl else []):
        command = [str(probe), str(eop), str(case/'input/orbit_keplerian.txt'), str(mode)]
        if mode == 1: command.append(str(args.jpl.resolve()))
        p = subprocess.run(command, capture_output=True, text=True, timeout=30)
        require(p.returncode==0, p.stderr)
        lines = p.stdout.splitlines()
        ae = float(lines[0])
        rows = [[float(x) for x in line.split()] for line in lines[1:]]
        require(len(rows)==8641 and all(len(r)==11 for r in rows), 'Malformed geometry probe')
        all_probe_rows[mode] = rows

        def run(name, changes=None, good=True):
            settings = {'object_id':'geometry', 'file_eph_Cartesian':'',
                        'file_eph_Keplerian':str(case/'input/orbit_keplerian.txt'),
                        'file_EOP':str(eop), 'is_input_normalized':'0',
                        'file_observer':str(case/'input/observers.txt'), 'observer_id':'1',
                        'need_output_Keplerian':'0', 'observation_type':'1',
                        'min_elevation_deg':'5', 'use_station_night':'0', 'use_earth_shadow':'0',
                        'max_sun_elevation_deg':'-6', 'solar_ephemeris':str(mode)}
            if mode == 1: settings['file_JPLEPH'] = str(args.jpl.resolve())
            settings.update(changes or {})
            out = root/name
            out.mkdir()
            settings['dir_output'] = str(out)
            config = root/(name+'.txt')
            config.write_text('\n'.join(k+' = '+v for k,v in settings.items())+'\n')
            p = subprocess.run([str(binary), '-I', str(config)],
                               capture_output=True, text=True, timeout=30)
            require((p.returncode == 0) == good, p.stderr[-3000:])

            if not good:
                require('[ERROR]' in p.stderr and 'Call Stack' in p.stderr, 'Noise failure lost tracer context')
                require('finished successfully' not in p.stderr, 'Noise failure reported success')
            return {float(line.split()[0]) for f in out.glob('obs_*.txt') for line in f.read_text().splitlines()}

        counts = []
        for night, shadow in [(0,0), (1,0), (0,1), (1,1)]:
            expected = set()
            for row in rows:
                _, utc = row[:2]
                station, target, sun = row[2:5], row[5:8], row[8:11]
                R, r = norm(station), norm(target)
                # 直接使用教材角度公式；生产代码的昼夜/地影判断采用等价点积式。
                high = angle(station, target) < math.acos(R/r*math.cos(math.radians(5)))-math.radians(5)
                dark = angle(station, sun) > math.pi/2+math.radians(6)
                lit = angle(target, sun) < math.pi-math.asin(1/r)
                if high and (not night or dark) and (not shadow or lit): expected.add(utc)
            actual = run(f'formula-{mode}-{night}-{shadow}', {'use_station_night':str(night), 'use_earth_shadow':str(shadow)})
            require(actual==expected, f'Textbook disagreement for mode/night/shadow {mode}/{night}/{shadow}: {len(actual)} vs {len(expected)}')
            counts.append(len(actual))

        # 完整文件中的无噪角度还须对不可见点成立，用站心矢量独立核对全部采样。
        name = f'all-theory-{mode}'
        run(name, {'need_output_all':'1'})
        all_file = next((root/name).glob('all_*.txt'))
        theory = [line.split() for line in all_file.read_text().splitlines()]
        require(len(theory)==len(rows), 'Complete theory does not cover every input record')
        largest_error = 0.0
        for result, reference in zip(theory, rows):
            dx, dy, dz = [a-b for a,b in zip(reference[5:8], reference[2:5])]
            ra = math.degrees(math.atan2(dy, dx)) % 360
            dec = math.degrees(math.atan2(dz, math.hypot(dx, dy)))
            largest_error = max(largest_error, abs((float(result[4])-ra+180)%360-180), abs(float(result[5])-dec))
        require(largest_error < 5e-10, f'Incorrect theoretical RA/Dec in complete file: {largest_error} deg')

        # 以探针中的实际测站和太阳方向构造边界两侧的几何，速度取0，关闭Kep真值。
        def single(name, row, position, changes=None):
            file = root/(name+'-input.txt')
            normalized = (changes or {}).get('is_input_normalized') == '1'
            xyz = position if normalized else [x*ae/1000 for x in position]
            file.write_text('1 '+format(row[1], '.17g')+' '+' '.join(format(x, '.17g') for x in xyz)+' 0 0 0\n')
            settings = {'file_eph_Keplerian':'', 'file_eph_Cartesian':str(file)}
            settings.update(changes or {})
            return len(run(name, settings))

        row = rows[0]
        station = row[2:5]
        n = unit(station)
        tangent = unit(cross(n, [0, 0, 1]))
        r = 2.0
        theta0 = math.acos(norm(station)/r*math.cos(math.radians(5)))-math.radians(5)
        for sign, expected in [(-1, 1), (1, 0)]:
            theta = theta0+sign*1e-8
            target = [r*(math.cos(theta)*a+math.sin(theta)*b) for a,b in zip(n,tangent)]
            require(single(f'height-{mode}-{sign}', row, target)==expected, 'Minimum-elevation boundary reversed')

        row = next(row for row in rows if -30 < math.degrees(math.asin(dot(unit(row[2:5]), unit(row[8:11])))) < -10)
        n, sun = unit(row[2:5]), unit(row[8:11])
        beta = math.degrees(math.asin(dot(n, sun)))
        target = [2*x for x in n]
        require(single(f'zenith-{mode}', row, target)==1, 'Overhead target was lost to acos roundoff')
        for offset, expected in [(1e-6, 1), (-1e-6, 0)]:
            require(single(f'night-{mode}-{offset}', row, target,
                           {'use_station_night':'1', 'max_sun_elevation_deg':format(beta+offset, '.17g')})==expected,
                    'Station-night boundary reversed')

        row = next(row for row in rows if dot(unit(row[2:5]), unit(row[8:11])) < -.5)
        sun = unit(row[8:11])
        tangent = unit(cross(sun, [0, 0, 1]))
        for offset, expected in [(-1e-7, 0), (1e-7, 1)]:
            target = [-10*a+(1+offset)*b for a,b in zip(sun,tangent)]
            require(single(f'shadow-{mode}-{offset}', row, target, {'use_earth_shadow':'1'})==expected,
                    'Cylindrical-shadow boundary reversed')
        require(single(f'normalized-cart-{mode}', row, target,
                       {'use_earth_shadow':'1', 'is_input_normalized':'1'})==1,
                'Normalized Cartesian input changed shadow visibility')
        normalized_kep = root/f'normalized-kep-{mode}-input.txt'
        data = []
        for line in (case/'input/orbit_keplerian.txt').read_text().splitlines():
            fields = line.split()
            elements = [float(x) for x in fields[2:]]
            elements[0] /= ae/1000
            elements[2:] = [math.radians(x) for x in elements[2:]]
            data.append(' '.join(fields[:2]+[format(x, '.17g') for x in elements]))
        normalized_kep.write_text('\n'.join(data)+'\n')
        require(len(run(f'normalized-kep-{mode}', {'file_eph_Keplerian':str(normalized_kep), 'is_input_normalized':'1'}))==counts[0],
                'Normalized Keplerian input changed visibility')
        summaries['JPL' if mode==1 else 'analytical'] = {'formula_counts_off_night_shadow_both':counts, 'boundary_cases':6, 'zenith_case':1, 'normalized_input_checks':2, 'all_theory_max_component_error_deg':largest_error}

    if len(all_probe_rows)==2:
        differences = [math.degrees(angle(a[8:11], b[8:11]))*3600
                       for a,b in zip(all_probe_rows[1], all_probe_rows[2])]
        summaries['maximum_Sun_direction_difference_arcsec'] = max(differences)
    print(json.dumps({'checks':checks, 'results':summaries}, ensure_ascii=False, indent=2))
