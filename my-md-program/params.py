"""アルゴンのLJパラメータとシミュレーション条件（すべてHartree原子単位系）。

SI単位系へ書き換えるときは、本ファイルの数値を kg / m / J / s に直し、
constants.py を SI 用に差し替えれば、potential.py 以下のロジックは変更不要。
"""

from constants import AMU_TO_AU_MASS, KB_AU

# ----- Ar 物性パラメータ -----
# 質量：39.948 amu ≈ 72820.7 m_e
MASS = 39.948 * AMU_TO_AU_MASS

# LJ σ：3.405 Å ≈ 6.4344 bohr
SIGMA = 6.4344

# LJ ε：119.8 K · kB ≈ 3.7938e-4 Hartree
EPSILON = 119.8 * KB_AU

# カットオフ距離：2.5σ
RC = 2.5 * SIGMA

# Cutoff energy shift：U(rc) を引いてエネルギーを連続にする（NVE保存のため必須）
_S_OVER_RC_6 = (SIGMA / RC) ** 6
_S_OVER_RC_12 = _S_OVER_RC_6 * _S_OVER_RC_6
E_SHIFT = 4.0 * EPSILON * (_S_OVER_RC_12 - _S_OVER_RC_6)

# ----- 計算条件 -----
N_CELLS = 4               # FCC 4×4×4 = 256 原子
LATTICE_A = 10.81         # bohr (≈5.72 Å, 液体Ar密度近傍)
T0 = 94.4                 # K（初期温度。Ar臨界点近傍の液相）
DT = 20.0                 # a.u. (≈0.48 fs)
N_STEPS = 1000
OUT_FREQ = 10
SEED = 42
