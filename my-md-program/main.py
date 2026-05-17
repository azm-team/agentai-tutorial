"""アルゴンLJ-MD（NVE）のエントリポイント。

実行：
    cd my-md-program
    python main.py

出力：
    trajectory.xyz  — VMD/Ovito等で開けるXYZ形式のトラジェクトリ
    energy.log      — step / time / T / KE / PE / E_tot / drift のログ
"""

from constants import BOHR_TO_ANGSTROM, AU_TIME_TO_FS, KB_AU
from params import (
    MASS, SIGMA, EPSILON, RC, E_SHIFT,
    N_CELLS, LATTICE_A, T0, DT, N_STEPS, OUT_FREQ, SEED,
)
from potential import lj_force_energy
from integrator import velocity_verlet_step
from initialize import init_fcc, init_velocities_mb
from thermo import kinetic_energy, temperature
from output import write_log_header, write_log_line, write_xyz_frame


def main():
    # 初期化
    pos, box_l = init_fcc(N_CELLS, LATTICE_A)
    N = pos.shape[0]
    vel = init_velocities_mb(N, MASS, T0, KB_AU, SEED)
    forces, pe = lj_force_energy(pos, box_l, SIGMA, EPSILON, RC, E_SHIFT)
    ke = kinetic_energy(vel, MASS)
    e0 = ke + pe

    print(
        f"原子数 N = {N},  "
        f"ボックス長 = {box_l:.3f} bohr ({box_l * BOHR_TO_ANGSTROM:.3f} Å)"
    )
    print(
        f"初期温度 T0 = {T0:.2f} K,  "
        f"dt = {DT:.2f} a.u. ({DT * AU_TIME_TO_FS:.4f} fs),  "
        f"steps = {N_STEPS}"
    )
    print(f"cutoff RC = {RC:.3f} bohr,  E_SHIFT = {E_SHIFT:.4e} Hartree")
    print(f"初期 E_tot = {e0:.6e} Hartree")
    print("-" * 60)

    # メインループ
    with open("trajectory.xyz", "w") as ftraj, open("energy.log", "w") as flog:
        write_log_header(flog)
        for step in range(N_STEPS + 1):
            if step > 0:
                pos, vel, forces, pe = velocity_verlet_step(
                    pos, vel, forces, box_l, MASS, DT,
                    SIGMA, EPSILON, RC, E_SHIFT,
                )

            if step % OUT_FREQ == 0:
                ke = kinetic_energy(vel, MASS)
                T = temperature(ke, N, KB_AU)
                e_tot = ke + pe
                drift = (e_tot - e0) / abs(e0)
                time_fs = step * DT * AU_TIME_TO_FS
                write_log_line(flog, step, time_fs, T, ke, pe, e_tot, drift)
                write_xyz_frame(ftraj, pos, step, BOHR_TO_ANGSTROM)
                if step % (OUT_FREQ * 10) == 0:
                    print(
                        f"step {step:5d}  T = {T:7.2f} K  "
                        f"E_tot = {e_tot:.6e}  drift = {drift:+.2e}"
                    )

    print("-" * 60)
    print("完了。trajectory.xyz と energy.log を出力しました。")


if __name__ == "__main__":
    main()
