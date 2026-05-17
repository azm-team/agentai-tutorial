"""Lennard-Jonesポテンシャルによる力とエネルギーの計算（@njitで高速化）。

物理量は単位系に依存しない（σ, ε, rc, e_shift を引数で受け取るため）。
立方体周期境界、最小イメージ規約、cutoff with energy shift を採用。
"""

import numpy as np
from numba import njit


@njit
def lj_force_energy(pos, box_l, sigma, epsilon, rc, e_shift):
    """LJ力と全ポテンシャルエネルギーを計算する。

    引数
        pos      : (N, 3) ndarray  原子座標
        box_l    : float           立方体ボックスの一辺
        sigma    : float           LJ σ
        epsilon  : float           LJ ε
        rc       : float           カットオフ距離（rc < box_l / 2 を仮定）
        e_shift  : float           cutoff energy shift（U(rc)）

    戻り値
        forces : (N, 3) ndarray  各原子に働く力
        pe     : float           全ポテンシャルエネルギー
    """
    N = pos.shape[0]
    forces = np.zeros_like(pos)
    pe = 0.0
    rc2 = rc * rc
    sig2 = sigma * sigma
    for i in range(N - 1):
        for j in range(i + 1, N):
            dx = pos[i, 0] - pos[j, 0]
            dy = pos[i, 1] - pos[j, 1]
            dz = pos[i, 2] - pos[j, 2]
            # 最小イメージ規約（立方体）
            dx -= box_l * np.rint(dx / box_l)
            dy -= box_l * np.rint(dy / box_l)
            dz -= box_l * np.rint(dz / box_l)
            r2 = dx * dx + dy * dy + dz * dz
            if r2 < rc2:
                inv_r2 = 1.0 / r2
                s2_r2 = sig2 * inv_r2
                s6 = s2_r2 * s2_r2 * s2_r2
                s12 = s6 * s6
                pe += 4.0 * epsilon * (s12 - s6) - e_shift
                # f_ij = -dU/dr · r̂ → ベクトル形：fmag_over_r * (r_i - r_j)
                fmag_over_r = 24.0 * epsilon * inv_r2 * (2.0 * s12 - s6)
                forces[i, 0] += fmag_over_r * dx
                forces[i, 1] += fmag_over_r * dy
                forces[i, 2] += fmag_over_r * dz
                forces[j, 0] -= fmag_over_r * dx
                forces[j, 1] -= fmag_over_r * dy
                forces[j, 2] -= fmag_over_r * dz
    return forces, pe
