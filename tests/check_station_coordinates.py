#!/usr/bin/env python3
"""核对六向测站文件转换、内存接口接续，以及输入错误的 Tracer。"""
import argparse
import json
import math
from pathlib import Path
import subprocess


def read_rows(path):
    """读取 UTC MJD 与三分量位置，忽略注释和空行。"""
    return [list(map(float, line.split('#', 1)[0].split()))
            for line in path.read_text().splitlines() if line.split('#', 1)[0].strip()]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--binary', type=Path, required=True)
    parser.add_argument('--probe', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    run = args.output.resolve()
    run.mkdir(parents=True, exist_ok=True)
    binary = args.binary.resolve()
    eop = root / 'data/eop/EOP-All_250516.txt'
    fixture = root / 'tests/convert_coordinates/geodetic.txt'

    # 科学参考由 C++ 内存接口产生；探针同时执行独立 SOFA 和极点检查。
    reference = run / 'memory_reference.txt'
    probe = subprocess.run([str(args.probe.resolve()), str(eop), str(reference)],
                           capture_output=True, text=True, timeout=30)
    (run / 'memory_probe.log').write_text(probe.stdout + probe.stderr)
    assert probe.returncode == 0, probe.stdout + probe.stderr
    memory = read_rows(reference)
    geodetic = read_rows(fixture)
    assert len(memory) == len(geodetic) == 9

    # 首两次转换准备地固与天球文件，随后验证另外四个方向。
    paths = {'geodetic': fixture, 'itrs': run / 'geodetic_to_itrs.txt',
             'gcrs': run / 'geodetic_to_gcrs.txt'}
    results = {}
    for source, target in [('geodetic', 'itrs'), ('geodetic', 'gcrs'),
                           ('itrs', 'geodetic'), ('itrs', 'gcrs'),
                           ('gcrs', 'itrs'), ('gcrs', 'geodetic')]:
        name = source + '_to_' + target
        output = run / (name + '.txt')
        command = [str(binary), source, target, str(paths[source]), str(output)]
        if 'gcrs' in (source, target):
            command.append(str(eop))
        process = subprocess.run(command, capture_output=True, text=True, timeout=30)
        (run / (name + '.log')).write_text(process.stdout + process.stderr)
        assert process.returncode == 0, (name, process.stderr)
        rows = read_rows(output)
        assert len(rows) == len(memory)
        assert [row[0] for row in rows] == [row[0] for row in geodetic], name
        assert all(math.isfinite(value) for row in rows for value in row), name

        # 直角坐标按分量对照内存输出；大地坐标检查高度、纬度及非极点经度。
        if target in ('itrs', 'gcrs'):
            offset = 1 if target == 'itrs' else 4
            difference = max(abs(row[axis + 1] - expected[axis + offset])
                             for row, expected in zip(rows, memory) for axis in range(3))
            assert difference < 1e-6, (name, difference)
            results[name] = {'max_component_difference_m': difference}
        else:
            height = max(abs(row[3] - expected[3]) for row, expected in zip(rows, geodetic))
            latitude = max(abs(row[2] - expected[2]) for row, expected in zip(rows, geodetic))
            longitude = max(abs(math.remainder(row[1] - expected[1], 360.0))
                            for row, expected in zip(rows, geodetic) if abs(expected[2]) < 89.999999)
            assert height < 1e-6 and latitude < 1e-10 and longitude < 1e-10, (name, height, latitude, longitude)
            results[name] = {'max_height_difference_m': height,
                             'max_latitude_difference_deg': latitude,
                             'max_longitude_difference_deg_excluding_poles': longitude}

    # 两种单位分别覆盖六向转换，换回米后与默认输入输出直接比较。
    unit_results = {}
    for units, factor in [('deg-m', 1.0), ('deg-km', 1000.0)]:
        unit_run = run / units
        unit_run.mkdir(exist_ok=True)
        inputs = {}
        for frame, path in paths.items():
            rows = read_rows(path)
            for row in rows:
                axes = [3] if frame == 'geodetic' else [1, 2, 3]
                for axis in axes:
                    row[axis] /= factor
            inputs[frame] = unit_run / (frame + '_input.txt')
            inputs[frame].write_text(''.join(' '.join(format(value, '.17g') for value in row) + '\n' for row in rows))

        unit_results[units] = {}
        for source, target in [('geodetic', 'itrs'), ('geodetic', 'gcrs'),
                               ('itrs', 'geodetic'), ('itrs', 'gcrs'),
                               ('gcrs', 'itrs'), ('gcrs', 'geodetic')]:
            name = source + '_to_' + target
            output = unit_run / (name + '.txt')
            command = [str(binary), source, target, str(inputs[source]), str(output)]
            if 'gcrs' in (source, target):
                command.append(str(eop))
            command += ['--units', units]
            process = subprocess.run(command, capture_output=True, text=True, timeout=30)
            (unit_run / (name + '.log')).write_text(process.stdout + process.stderr)
            assert process.returncode == 0, (units, name, process.stderr)
            rows = read_rows(output)
            expected = read_rows(run / (name + '.txt'))
            assert len(rows) == len(expected)
            assert [row[0] for row in rows] == [row[0] for row in expected]
            axes = [3] if target == 'geodetic' else [1, 2, 3]
            difference = max(abs(row[axis] * factor - ref[axis])
                             for row, ref in zip(rows, expected) for axis in axes)
            assert difference < 1e-6, (units, name, difference)
            if target == 'geodetic':
                assert all(abs(row[2] - ref[2]) < 1e-10 for row, ref in zip(rows, expected))
                assert all(abs(math.remainder(row[1] - ref[1], 360.0)) < 1e-10
                           for row, ref in zip(rows, expected) if abs(ref[2]) < 89.999999)
            assert ('height_' if target == 'geodetic' else 'X_') + ('m' if factor == 1.0 else 'km') in output.read_text().splitlines()[0]
            if units == 'deg-m':
                assert output.read_bytes() == (run / (name + '.txt')).read_bytes(), name
            unit_results[units][name] = {'max_length_component_difference_m': difference}

    # 所有输入均在写出前检查；坏输入应保留已有结果，错误中应有主程序调用链。
    protected = run / 'preserved.txt'
    protected.write_text('preserved output\n')
    malformed = {
        'extra_column': '60000 1 2 3 4\n',
        'missing_column': '60000 1 2\n',
        'nonfinite': '60000 1 nan 3\n',
        'empty': '# no records\n',
        'latitude': '60000 10 91 0\n',
        'late_bad_record': '60000 10 20 0\n60000 10 91 0\n'
    }
    bad_commands = []
    for name, text in malformed.items():
        path = run / (name + '.txt')
        path.write_text(text)
        bad_commands.append((name, ['geodetic', 'itrs', str(path), str(protected)]))
    origin = run / 'origin.txt'
    origin.write_text('60000 0 0 0\n')
    beyond = run / 'outside_eop.txt'
    beyond.write_text('90000 0 0 0\n')
    bad_commands += [
        ('unsupported_units', ['geodetic', 'itrs', str(fixture), str(protected), '--units', 'rad-m']),
        ('unsupported_length', ['geodetic', 'itrs', str(fixture), str(protected), '--units', 'deg-cm']),
        ('missing_units_value', ['geodetic', 'itrs', str(fixture), str(protected), '--units']),
        ('zero_itrs', ['itrs', 'geodetic', str(origin), str(protected)]),
        ('same_frame', ['itrs', 'itrs', str(paths['itrs']), str(protected)]),
        ('unknown_frame', ['unknown', 'itrs', str(fixture), str(protected)]),
        ('missing_eop', ['geodetic', 'gcrs', str(fixture), str(protected)]),
        ('extra_eop', ['geodetic', 'itrs', str(fixture), str(protected), str(eop)]),
        ('bad_eop_path', ['geodetic', 'gcrs', str(fixture), str(protected), str(run / 'missing_eop.txt')]),
        ('outside_eop', ['itrs', 'gcrs', str(beyond), str(protected), str(eop)]),
        ('missing_input', ['geodetic', 'itrs', str(run / 'missing.txt'), str(protected)]),
        ('missing_output_parent', ['geodetic', 'itrs', str(fixture), str(run / 'missing_dir/output.txt')])
    ]
    if Path('/dev/full').exists():
        bad_commands.append(('output_full', ['geodetic', 'itrs', str(fixture), '/dev/full']))
    for name, command in bad_commands:
        process = subprocess.run([str(binary)] + command, capture_output=True, text=True, timeout=10)
        (run / (name + '.log')).write_text(process.stdout + process.stderr)
        assert process.returncode != 0, name
        assert 'main_convert_coordinates.cpp' in process.stderr, (name, process.stderr)
        assert protected.read_text() == 'preserved output\n', name

    # 版本号来自整库 CMake，帮助与批量入口均无需交互。
    for option in ('--help', '--version'):
        process = subprocess.run([str(binary), option], capture_output=True, text=True, timeout=10)
        assert process.returncode == 0 and 'aster_convert_coordinates' in process.stdout
        (run / (option[2:] + '.txt')).write_text(process.stdout)
    summary = {'baseline': 'C++ memory API; WGS84 and SOFA IAU2000B assembly checked by probe',
               'records': len(memory), 'directions': results, 'units': unit_results, 'failure_cases': len(bad_commands),
               'iod_three_observation_gcrs_m': read_rows(paths['gcrs'])[:3]}
    (run / 'summary.json').write_text(json.dumps(summary, ensure_ascii=False, indent=2) + '\n')
    print('PASS six file directions, IOD memory handoff, round trips and', len(bad_commands), 'failure cases')


if __name__ == '__main__':
    main()
