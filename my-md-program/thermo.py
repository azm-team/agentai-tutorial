"""熱力学量（運動エネルギー、温度、重心速度）の計算。"""

import numpy as np


def kinetic_energy(vel, mass):
    """全運動エネルギー KE = 0.5 * m * Σ|v|²（全原子同質量を仮定）。"""
    return 0.5 * mass * np.sum(vel * vel)


def temperature(ke, N, kB):
    """瞬時温度。重心並進3自由度を除いた dof = 3N - 3 で算出。"""
    return 2.0 * ke / ((3 * N - 3) * kB)


def com_velocity(vel):
    """重心速度（COM drift の診断用）。"""
    return vel.mean(axis=0)


def remove_com_drift(vel):
    """重心並進運動を除去する。"""
    vel -= vel.mean(axis=0)
    return vel
