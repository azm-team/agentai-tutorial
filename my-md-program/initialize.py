"""初期配置（FCC格子）と初期速度（Maxwell-Boltzmann）の生成。"""

import numpy as np


def init_fcc(n_cells, a):
    """FCC格子に原子を配置する。

    1ユニットセル4原子の basis：
        (0, 0, 0), (0.5, 0.5, 0), (0.5, 0, 0.5), (0, 0.5, 0.5)  （× a）

    引数
        n_cells : int    1辺あたりのユニットセル数
        a       : float  ユニットセルの格子定数

    戻り値
        pos   : (4*n_cells**3, 3) ndarray  原子座標
        box_l : float                       シミュレーションボックスの一辺（= n_cells * a）
    """
    basis = np.array([
        [0.0, 0.0, 0.0],
        [0.5, 0.5, 0.0],
        [0.5, 0.0, 0.5],
        [0.0, 0.5, 0.5],
    ])
    N = 4 * n_cells ** 3
    pos = np.empty((N, 3))
    idx = 0
    for ix in range(n_cells):
        for iy in range(n_cells):
            for iz in range(n_cells):
                cell_origin = np.array([ix, iy, iz], dtype=float)
                for b in basis:
                    pos[idx] = (cell_origin + b) * a
                    idx += 1
    box_l = n_cells * a
    return pos, box_l


def init_velocities_mb(N, mass, T, kB, seed):
    """Maxwell-Boltzmann分布で速度を生成し、COMドリフト除去後に温度 T に厳密リスケール。

    引数
        N    : int    原子数
        mass : float  1原子あたりの質量
        T    : float  目標温度
        kB   : float  Boltzmann定数
        seed : int    乱数seed

    戻り値
        vel : (N, 3) ndarray  初期速度
    """
    rng = np.random.default_rng(seed)
    sigma_v = np.sqrt(kB * T / mass)
    vel = rng.normal(0.0, sigma_v, size=(N, 3))

    # COMドリフトを除去（全体並進運動量をゼロに）
    vel -= vel.mean(axis=0)

    # 残った自由度（3N - 3）で温度を計算して目標温度に厳密リスケール
    ke = 0.5 * mass * np.sum(vel * vel)
    dof = 3 * N - 3
    T_current = 2.0 * ke / (dof * kB)
    vel *= np.sqrt(T / T_current)
    return vel
