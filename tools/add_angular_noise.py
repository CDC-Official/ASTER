#!/usr/bin/env python3
"""给理论赤经、赤纬 OBS 加入独立的东西、南北高斯误差，并生成蒙特卡洛样本。

输入为 C++ 生成器写出的六列 OBS：UTC MJD、类型、目标号、站号、RA(deg)、Dec(deg)。
每组样本沿用输入的时刻、身份和分弧；内部使用 rad，输出保持 deg。
仅使用 Python 标准库。示例：
    python3 tools/add_angular_noise.py theory/obs_*.txt --sigma-arcsec 1 --samples 100 --seed 42 --output noisy
"""

import argparse
import hashlib
import json
import math
from pathlib import Path
import platform
import random
import shutil
import sys

# 与工程 ARCSEC_TO_RAD 的定义一致；标准差输入为 arcsec，角度计算使用 rad。
ARCSEC_TO_RAD = math.pi / (180.0 * 3600.0)


def read_theory(path):
    """读取一条理论弧，返回原始六列字符串及 rad 制赤经、赤纬。"""
    records = []
    previous_time = -math.inf

    for line_number, line in enumerate(path.read_text().splitlines(), 1):
        if not line.strip() or line.lstrip().startswith('#'):
            continue

        fields = line.split()

        if len(fields) != 6 or fields[1] != '1':
            raise ValueError(f'{path}:{line_number}: expected six-column ground RA/Dec OBS (type 1)')

        utc, ra_deg, dec_deg = (float(fields[k]) for k in (0, 4, 5))

        if not all(math.isfinite(x) for x in (utc, ra_deg, dec_deg)):
            raise ValueError(f'{path}:{line_number}: non-finite epoch or angle')

        if utc < previous_time or not 0 <= ra_deg < 360 or not -90 <= dec_deg <= 90:
            raise ValueError(f'{path}:{line_number}: invalid epoch order or angular domain')

        records.append((fields, math.radians(ra_deg), math.radians(dec_deg)))
        previous_time = utc

    if not records:
        raise ValueError(f'{path}: empty theoretical arc')

    return records


def add_noise(records, sigma_rad, generator):
    """在完整弧上逐点抽取误差，返回六列含噪 OBS 文本；两个方向共用 sigma_rad。"""
    lines = []

    for record_number, (fields, ra0, dec0) in enumerate(records, 1):
        # 赤经改变量使用理论赤纬；赤道坐标极点处该小角度模型没有定义。
        if abs(dec0) >= math.pi / 2:
            raise ValueError(f'record {record_number}: east/north small-angle noise is undefined at a coordinate pole')

        east_error = generator.gauss(0.0, sigma_rad)
        north_error = generator.gauss(0.0, sigma_rad)
        ra = ra0 + east_error / math.cos(dec0)
        dec = dec0 + north_error

        # 不截断分母，不反射赤纬，不重抽误差；超出模型定义域时报告本条记录。
        if not math.isfinite(ra) or not math.isfinite(dec) or abs(dec) >= math.pi / 2:
            raise ValueError(f'record {record_number}: non-finite angle or noisy declination reached a pole')

        ra_deg = math.degrees(ra % math.tau)

        if ra_deg >= 360:
            ra_deg = 0.0

        lines.append(' '.join([*fields[:4], format(ra_deg, '.17g'), format(math.degrees(dec), '.17g')]))

    return '\n'.join(lines) + '\n'


def main():
    """读入设置与理论弧，逐组写出含噪 OBS，并保存随机数及源资料说明。"""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('inputs', type=Path, nargs='+', help='one or more theoretical six-column OBS files')
    parser.add_argument('--sigma-arcsec', type=float, required=True, help='common east/north Gaussian standard deviation')
    parser.add_argument('--samples', type=int, default=1, help='number of independent realizations, default 1')
    parser.add_argument('--seed', type=int, required=True, help='nonnegative base random seed')
    parser.add_argument('--output', type=Path, required=True, help='new output directory')
    args = parser.parse_args()

    try:
        # 检查全部输入后再创建结果目录；同名输入会使每组的输出文件冲突。
        sigma_rad = args.sigma_arcsec * ARCSEC_TO_RAD

        if not math.isfinite(args.sigma_arcsec) or args.sigma_arcsec < 0 or not math.isfinite(sigma_rad):
            raise ValueError('--sigma-arcsec must be finite and nonnegative')

        if args.sigma_arcsec > 0 and sigma_rad == 0:
            raise ValueError('--sigma-arcsec is too small to represent in radians')

        if args.samples < 1 or args.seed < 0:
            raise ValueError('--samples must be positive and --seed must be nonnegative')

        sources = [path.resolve(strict=True) for path in args.inputs]

        if len({path.name for path in sources}) != len(sources):
            raise ValueError('input filenames must be distinct')

        arcs = []

        for path in sources:
            records = read_theory(path)
            digest = hashlib.sha256(path.read_bytes()).hexdigest()

            if sigma_rad > 0 and any(abs(dec0) >= math.pi / 2 for _, _, dec0 in records):
                raise ValueError(f'{path}: positive noise is undefined at a coordinate pole')

            arcs.append((path, records, digest))

        # 源资料只读，每组保留原弧文件名；轨道真值继续保存在理论资料目录。
        args.output.mkdir(parents=True, exist_ok=False)
        manifest = {
            'model': 'independent_equal_east_north_small_angle',
            'formula': 'RA=wrap(RA0+eta_E/cos(Dec0)); Dec=Dec0+eta_N',
            'sigma_arcsec': args.sigma_arcsec,
            'samples': args.samples,
            'base_seed': args.seed,
            'generator': 'Python random.Random / MT19937 / gauss',
            'python_version': platform.python_version(),
            'seed_derivation': 'SHA256 of JSON [base_seed, source_sha256, sample_number]; sample_number starts at 1',
            'script_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            'sources': [{'path': str(path), 'sha256': digest, 'records': len(records)} for path, records, digest in arcs],
        }
        (args.output / 'noise.json').write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + '\n')

        for sample_number in range(1, args.samples + 1):
            sample_dir = args.output / f'sample_{sample_number:06d}'
            sample_dir.mkdir()

            for path, records, digest in arcs:
                output_file = sample_dir / path.name

                # 零噪声直接复制，保留原始精度、空行和注释。
                if sigma_rad == 0:
                    shutil.copyfile(path, output_file)
                    continue

                # 用内容摘要稳定识别输入弧；移动文件、调整输入顺序或增加样本组数不改变已有结果。
                seed_text = json.dumps([args.seed, digest, sample_number], separators=(',', ':'))
                derived_seed = int.from_bytes(hashlib.sha256(seed_text.encode()).digest(), 'big')
                generator = random.Random(derived_seed)

                try:
                    text = add_noise(records, sigma_rad, generator)
                except ValueError as error:
                    raise ValueError(f'{path}, sample {sample_number}: {error}') from error

                output_file.write_text(text)

        print(f'Generated {args.samples} samples of {len(arcs)} arcs in {args.output.resolve()}')
        return 0

    except (OSError, ValueError, OverflowError) as error:
        parser.exit(1, f'Error: {error}\n')


if __name__ == '__main__':
    sys.exit(main())
