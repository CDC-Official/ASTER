#!/usr/bin/env python3
"""小规模核对 Python 测角加噪的定义、复现性、身份列和坐标边界。"""
import argparse
import hashlib
import json
import math
from pathlib import Path
import random
import statistics
import subprocess
import sys
import tempfile


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--script', type=Path, required=True)
    args = parser.parse_args()
    script = args.script.resolve()
    sigma_rad = math.pi / (180 * 3600)
    summary = {}

    with tempfile.TemporaryDirectory(prefix='cdc-angular-noise-') as folder:
        root = Path(folder)
        sources = []

        # 完整弧中使用普通赤纬和南北近极区，赤经置于 360 度左侧以覆盖跨界。
        for label, dec in [('ordinary', 35), ('north', 89.9), ('south', -89.9)]:
            source = root / (label + '.obs')
            source.write_text(''.join(f'{60000+k/86400:.17g} 1 TARGET STATION 359.999999 {dec}\n' for k in range(4096)))
            sources.append(source)

        def run(name, paths=sources, sigma=1, samples=2, seed=42, good=True):
            output = root / name
            command = [sys.executable, str(script), *map(str, paths), '--sigma-arcsec', str(sigma),
                       '--samples', str(samples), '--seed', str(seed), '--output', str(output)]
            result = subprocess.run(command, capture_output=True, text=True, timeout=20)
            assert (result.returncode == 0) == good, result.stdout + result.stderr
            return output

        noisy = run('noisy')
        repeated = run('repeat', paths=list(reversed(sources)), samples=3)
        zero = run('zero', sigma=0, samples=1)
        changed = run('changed', seed=43, samples=1)
        manifest = json.loads((noisy / 'noise.json').read_text())
        assert manifest['sigma_arcsec'] == 1 and len(manifest['sources']) == 3

        for source in sources:
            reference = [line.split() for line in source.read_text().splitlines()]
            output = noisy / 'sample_000001' / source.name
            samples = [line.split() for line in output.read_text().splitlines()]
            assert output.read_bytes() == (repeated / 'sample_000001' / source.name).read_bytes()
            assert (noisy / 'sample_000002' / source.name).read_bytes() == (repeated / 'sample_000002' / source.name).read_bytes()
            assert output.read_bytes() != (noisy / 'sample_000002' / source.name).read_bytes()
            assert output.read_bytes() != (changed / 'sample_000001' / source.name).read_bytes()
            assert source.read_bytes() == (zero / 'sample_000001' / source.name).read_bytes()
            assert len(samples) == len(reference) and all(a[:4] == b[:4] for a, b in zip(samples, reference))

            # 从写出的角度独立恢复两个方向的误差，不调用待测脚本的计算函数。
            east, north = [], []

            for observed, theory in zip(samples, reference):
                ra, dec, ra0, dec0 = map(math.radians, map(float, [*observed[4:], *theory[4:]]))
                east.append(math.remainder(ra-ra0, math.tau)*math.cos(dec0)/sigma_rad)
                north.append((dec-dec0)/sigma_rad)
                assert 0 <= float(observed[4]) < 360 and -90 < float(observed[5]) < 90

            means = [statistics.mean(east), statistics.mean(north)]
            deviations = [statistics.stdev(east), statistics.stdev(north)]
            assert max(map(abs, means)) < 5 / math.sqrt(len(samples)), means
            assert max(abs(x-1) for x in deviations) < 5 / math.sqrt(2*(len(samples)-1)), deviations
            correlation = sum((e-means[0])*(n-means[1]) for e, n in zip(east, north)) / ((len(samples)-1)*deviations[0]*deviations[1])
            assert abs(correlation) < 5 / math.sqrt(len(samples))
            assert any(float(row[4]) < 1 for row in samples) and any(float(row[4]) > 359 for row in samples)

            # 独立重建首条抽样，精确检查单位、抽样顺序和理论赤纬分母。
            digest = hashlib.sha256(source.read_bytes()).hexdigest()
            seed_text = json.dumps([42, digest, 1], separators=(',', ':'))
            seed = int.from_bytes(hashlib.sha256(seed_text.encode()).digest(), 'big')
            generator = random.Random(seed)
            assert abs(east[0]-generator.gauss(0, 1)) < 5e-8
            assert abs(north[0]-generator.gauss(0, 1)) < 5e-8
            summary[source.stem] = {'mean_arcsec': means, 'sd_arcsec': deviations, 'samples': len(samples)}

        # 极点允许零噪声原样复制；正噪声和越极点明确失败，不进行反射或重抽。
        pole = root / 'pole.obs'
        pole.write_text('60000 1 T S 0 90\n')
        run('pole-zero', [pole], sigma=0, samples=1)
        run('pole-noisy', [pole], good=False)
        crossing = root / 'crossing.obs'
        crossing.write_text(''.join(f'{60000+k/86400:.17g} 1 T S 0 89.999999\n' for k in range(64)))
        run('crossing-noisy', [crossing], good=False)
        run('invalid-sigma', sigma=-1, good=False)
        run('invalid-count', samples=0, good=False)
        run('invalid-seed', seed=-1, good=False)
        run('noisy', good=False)

    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == '__main__':
    main()
