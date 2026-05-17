# my-md-program

アルゴン原子の分子動力学（MD）シミュレーションを行う **教育用** Pythonコード。
Lennard-Jonesポテンシャル + Velocity Verlet + NVEアンサンブル。
高速化に [numba](https://numba.pydata.org/) の `@njit` を使用。

## 特徴

- 純粋な numpy + numba 実装（ASE等の高レベルライブラリ不使用）
- **Hartree原子単位系**（bohr, hartree, m_e, ℏ/Eh）で記述
- 周期境界（立方体）+ 最小イメージ規約 + cutoff energy shift
- 数秒〜十数秒で 1000 ステップ完走（256原子）

## 実行方法

```bash
# 依存：numpy, numba
pip install numpy numba

cd my-md-program
python main.py
```

出力：

- `trajectory.xyz` — XYZ形式のトラジェクトリ（座標はÅに変換）。VMD/Ovito等で可視化可能
- `energy.log` — `step / time[fs] / T[K] / KE / PE / E_tot / (E-E0)/|E0|`

## ファイル構成

| ファイル | 役割 |
|---|---|
| `main.py` | エントリポイント。初期化 → メインループ → 出力 |
| `constants.py` | **単位変換のみ**（bohr↔Å, Hartree↔eV, a.u.時間↔fs, kB） |
| `params.py` | Ar物性パラメータ（σ, ε, mass）と計算条件（dt, T0, N_CELLS, …） |
| `potential.py` | LJ力・ポテンシャルエネルギー（`@njit`） |
| `integrator.py` | Velocity Verlet 1ステップ（`@njit`） |
| `initialize.py` | FCC格子初期配置、Maxwell-Boltzmann速度 |
| `thermo.py` | 運動エネルギー、温度、重心速度 |
| `output.py` | XYZフレームとエネルギーログの書き出し |

## 設計メモ

**単位系と物理計算の分離**：

- `constants.py` は単位変換だけ、`params.py` は物質パラメータだけを持つ
- `potential.py` 以下のロジックは σ/ε/質量を**引数で受け取る**ため、単位系に依存しない
- → SI単位系へ書き換えるときは `constants.py` を SI 用に差し替え、`params.py` の数値を kg/m/J に直すだけで済む（`potential.py`, `integrator.py` 等は触らなくてよい）

## 動作確認

```bash
python main.py
# energy.log の (E-E0)/|E0| 列が ±1e-3 以内に収まっていれば NVE エネルギー保存 OK
python -c "import numpy as np; d=np.loadtxt('energy.log'); print('max |drift| =', np.max(np.abs(d[:,6])))"
```

NVE では運動エネルギーとポテンシャルエネルギーがやり取りされ、全エネルギーは
小さな数値誤差以外で保存する。初期温度 94.4 K の Maxwell-Boltzmann 分布で開始したあと、
equipartition により ~250 ステップで T ≈ 47 K 付近で揺らぐようになる
（初期PEがほぼ0で、平衡で半分がPEに移るため）。
