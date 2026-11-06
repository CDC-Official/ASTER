#!/usr/bin/env python3
"""用相同物理初态和独立DOP853积分，比较两套原生ASTER的FG预报。"""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess

import numpy as np
import scipy
from scipy.integrate import solve_ivp


def integrate(state, times, j2, period, tighter):
    """对归一化二体加J2方程正反向积分，保留所有输入时刻。"""
    def equation(_, current):
        position = current[:3]
        radius2 = position @ position
        radius = np.sqrt(radius2)
        fraction = position[2] ** 2 / radius2
        acceleration = -position / radius ** 3
        acceleration += 1.5 * j2 / radius ** 5 * position * np.array([
            5 * fraction - 1, 5 * fraction - 1, 5 * fraction - 3])
        return np.r_[current[3:], acceleration]

    result = np.empty((len(times), 3))

    for positive in (False, True):
        mask = times >= 0 if positive else times < 0
        selected = times[mask] if positive else times[mask][::-1]
        solution = solve_ivp(equation, (0, selected[-1]), state, method="DOP853", t_eval=selected,
                             rtol=5e-14 if tighter else 3e-13,
                             atol=5e-16 if tighter else 3e-15,
                             max_step=period / (512 if tighter else 256))

        if not solution.success:
            raise RuntimeError(solution.message)

        result[mask] = solution.y[:3].T if positive else solution.y[:3].T[::-1]

    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--before", type=Path, required=True)
    parser.add_argument("--after", type=Path, required=True)
    parser.add_argument("--case", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True, help="新的结果目录。")
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    record = json.loads(args.case.read_text())
    constants = record["constants"]
    case = record["case"]
    state = np.array(case["state_norm"])
    length = constants["unit_length_m"]
    time_scale = constants["unit_time_s"]
    scales = np.array([length] * 3 + [length / time_scale] * 3)
    times = np.arange(-300.0, 301.0)
    native_input = args.output / "prediction.txt"
    native_input.write_text("1\n" + f"{case['case_id']} 2 60400 {len(times)} 0.1\n"
                            + " ".join(f"{x:.17g}" for x in state * scales) + "\n"
                            + " ".join(f"{t:.17g}" for t in times) + "\n")
    predictions = []
    builds = []

    for role, binary in (("original", args.before), ("corrected", args.after)):
        # 先核对两套程序的物理常数，避免把参数差异混入系数修正。
        binary = binary.resolve()
        actual = json.loads(subprocess.check_output([str(binary), "--constants"], text=True))
        assert actual == constants, (role, actual, constants)
        path = args.output / (role + ".csv")
        subprocess.run([str(binary), str(native_input.resolve()), str(path.resolve())], check=True)
        predictions.append(np.loadtxt(path, delimiter=",", skiprows=1, usecols=(2, 3, 4)))
        builds.append({"role": role, "binary": str(binary),
                       "sha256": hashlib.sha256(binary.read_bytes()).hexdigest()})

    reference = integrate(state, times / time_scale, constants["j2_normalized"], case["period_s"] / time_scale, True) * length
    ordinary = integrate(state, times / time_scale, constants["j2_normalized"], case["period_s"] / time_scale, False) * length
    errors = [np.linalg.norm(prediction - reference, axis=1) for prediction in predictions]
    refinement = np.linalg.norm(ordinary - reference, axis=1)
    np.savetxt(args.output / "prediction.csv", np.column_stack([times, *errors, refinement]),
               delimiter=",", comments="", fmt="%.17g",
               header="time_s,original_position_difference_m,corrected_position_difference_m,reference_refinement_difference_m")
    summary = {"case": case, "constants": constants, "executables": builds,
               "scipy": scipy.__version__, "numpy": np.__version__,
               "case_sha256": hashlib.sha256(args.case.read_bytes()).hexdigest(),
               "max_reference_refinement_m": float(refinement.max()),
               "at_300_s_m": {"original": float(errors[0][-1]), "corrected": float(errors[1][-1])}}
    (args.output / "summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2) + "\n")
    print(json.dumps(summary["at_300_s_m"], ensure_ascii=False))


if __name__ == "__main__":
    main()
