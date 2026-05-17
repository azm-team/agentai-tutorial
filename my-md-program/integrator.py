"""Velocity Verlet時間積分（@njitで高速化）。

1ステップを1関数（kick → drift → 力再計算 → kick）で完結させることで、
順序ミスを防ぎ、AIに読ませた際の説明も安定する。
"""

import numpy as np
from numba import njit

from potential import lj_force_energy


@njit
def velocity_verlet_step(pos, vel, forces, box_l, mass, dt,
                          sigma, epsilon, rc, e_shift):
    """Velocity Verletで1ステップ進める。

    手順：half-kick → drift（PBCラップ込み）→ 力の再計算 → half-kick

    戻り値
        pos, vel, forces : 更新後の配列（in-placeでも更新される）
        pe               : 更新後のポテンシャルエネルギー
    """
    half_dt_over_m = 0.5 * dt / mass

    # half-kick
    vel += half_dt_over_m * forces

    # drift（PBCで [0, box_l) にラップ）
    pos += dt * vel
    pos -= box_l * np.floor(pos / box_l)

    # 力の再計算
    forces, pe = lj_force_energy(pos, box_l, sigma, epsilon, rc, e_shift)

    # half-kick
    vel += half_dt_over_m * forces

    return pos, vel, forces, pe
