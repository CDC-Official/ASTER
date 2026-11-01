#!/usr/bin/env python3
"""核对初轨文件与显式初始状态进入批处理后的结果一致。"""
import argparse
from pathlib import Path
import subprocess
from run_case import prepare


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--binary', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    tests = Path(__file__).resolve().parent
    for kind in ['gl', 'glml']:
        solution = tests / ('iod_' + kind) / 'reference' / (kind + '_solution.txt')
        fields = dict(line.split(' = ', 1) for line in solution.read_text().splitlines() if ' = ' in line)
        residuals = []
        for mode in ['file', 'inline']:
            _, _, run, argv = prepare('odn_batch', args.output.resolve() / (kind + '-' + mode))
            setting = run / 'setting.txt'
            values = dict(line.split(' = ', 1) for line in setting.read_text().splitlines() if ' = ' in line)
            values.update(file_observation=str(tests / 'iod_gl/input/observations.obs'),
                          file_observer=str(tests / 'iod_gl/input/observers.txt'),
                          observations_to_use='30', noise_sigma_arcsecond='3', iteration_maxround='2',
                          apriori_Keplerian='', file_orbit_solution_BODN='')
            if mode == 'file':
                values.update(file_apriori_orbit=str(solution), apriori_epoch='', apriori_Cartesian='')
            else:
                values.update(file_apriori_orbit='', apriori_epoch=fields['reference_epoch'],
                              apriori_Cartesian=fields['reference_Cartesian'])
            setting.write_text(''.join(key + ' = ' + value + '\n' for key, value in values.items()))
            proc = subprocess.run([str(args.binary.resolve())] + argv, cwd=run, capture_output=True, timeout=45)
            (run / 'console.txt').write_bytes(proc.stdout + proc.stderr)
            assert proc.returncode in [0, 2], (kind, mode, proc.stderr)
            residuals.append((run / 'bodn_residuals.txt').read_bytes())
        assert residuals[0] == residuals[1], kind
    print('PASS: GL and GLML file/inline initial states give identical batch residuals')


if __name__ == '__main__':
    main()
