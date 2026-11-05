#!/usr/bin/env python3
"""仅检查位置坐标和轨道状态在度米、度千米两种单位下的往返转换。"""
import argparse
import json
import math
from pathlib import Path
import subprocess


def load(path):
    """读取空白分列的数值记录，忽略 # 注释。"""
    return [list(map(float, line.split('#', 1)[0].split()))
            for line in path.read_text().splitlines() if line.split('#', 1)[0].strip()]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--bin-dir', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    binary = args.bin_dir.resolve()
    summary = {}

    for units, factor in [('deg-m', 1.0), ('deg-km', 1000.0)]:
        run = output / units
        run.mkdir(exist_ok=True)
        serial = 0

        def convert(program, source, target, input_path):
            """调用当前入口，保存一条转换命令的结果和日志。"""
            nonlocal serial
            serial += 1
            result = run / f'{serial:02d}_{source}_to_{target}.txt'
            command = [str(binary / program), source, target, str(input_path), str(result)]
            if program == 'aster_convert_coordinates' and 'gcrs' in (source, target):
                command.append(str(root / 'data/eop/EOP-All_250516.txt'))
            command += ['--units', units]
            process = subprocess.run(command, capture_output=True, text=True, timeout=15)
            result.with_suffix('.log').write_text(process.stdout + process.stderr)
            assert process.returncode == 0, process.stderr
            rows = load(result)
            assert rows and all(math.isfinite(value) for row in rows for value in row)
            assert [row[0] for row in rows] == [row[0] for row in load(input_path)]
            return result

        # 三对位置坐标各往返一次，包含原算例中的极点与近极点。
        coordinates = 'aster_convert_coordinates'
        geodetic = root / ('tests/convert_coordinates/geodetic.txt' if factor == 1.0
                           else 'tests/convert_coordinates/geodetic_km.txt')
        max_position_m = 0.0
        for intermediate in ['itrs', 'gcrs']:
            converted = convert(coordinates, 'geodetic', intermediate, geodetic)
            recovered = convert(coordinates, intermediate, 'geodetic', converted)
            for original, row in zip(load(geodetic), load(recovered)):
                max_position_m = max(max_position_m, abs(original[3] - row[3]) * factor)
                assert abs(original[2] - row[2]) < 1e-9
                if abs(original[2]) < 89.999999:
                    assert abs(math.remainder(original[1] - row[1], 360.0)) < 1e-9
        itrs = run / '01_geodetic_to_itrs.txt'
        gcrs = convert(coordinates, 'itrs', 'gcrs', itrs)
        recovered = convert(coordinates, 'gcrs', 'itrs', gcrs)
        max_position_m = max(max_position_m, max(abs(a - b) * factor
                             for first, second in zip(load(itrs), load(recovered))
                             for a, b in zip(first[1:], second[1:])))
        assert max_position_m < 1e-6, (units, max_position_m)

        # 两条普通椭圆轨道覆盖 Cartesian、Keplerian、NSG2 三对表示的往返。
        elements = load(root / 'tests/convert_states/keplerian_km.txt')
        for row in elements:
            row[1] *= 1000.0 / factor
        kep = run / 'keplerian_input.txt'
        kep.write_text(''.join(' '.join(format(value, '.17g') for value in row) + '\n' for row in elements))
        states = 'aster_convert_states'
        max_a_m = 0.0
        cartesian = None
        for intermediate in ['cart', 'nsg2']:
            converted = convert(states, 'kep', intermediate, kep)
            if intermediate == 'cart':
                cartesian = converted
            recovered = convert(states, intermediate, 'kep', converted)
            for original, row in zip(elements, load(recovered)):
                max_a_m = max(max_a_m, abs(original[1] - row[1]) * factor)
                assert abs(original[2] - row[2]) < 1e-11
                assert all(abs(math.remainder(a - b, 360.0)) < 1e-8 for a, b in zip(original[3:], row[3:]))
        nsg2 = convert(states, 'cart', 'nsg2', cartesian)
        recovered = convert(states, 'nsg2', 'cart', nsg2)
        pairs = list(zip(load(cartesian), load(recovered)))
        max_cart_m = max(abs(a - b) * factor for first, second in pairs for a, b in zip(first[1:4], second[1:4]))
        max_velocity_m_s = max(abs(a - b) * factor for first, second in pairs for a, b in zip(first[4:], second[4:]))
        assert max_a_m < 1e-5 and max_cart_m < 1e-5 and max_velocity_m_s < 1e-8
        summary[units] = {'coordinate_round_trip_max_length_difference_m': max_position_m,
                          'state_round_trip_max_semimajor_axis_difference_m': max_a_m,
                          'state_round_trip_max_position_component_difference_m': max_cart_m,
                          'state_round_trip_max_velocity_component_difference_m_s': max_velocity_m_s}

    (output / 'summary.json').write_text(json.dumps(summary, indent=2) + '\n')
    print('PASS coordinate and state round trips in deg-m and deg-km')


if __name__ == '__main__':
    main()
