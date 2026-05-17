"""トラジェクトリ（XYZ形式）とエネルギーログの出力。

ファイル名を ``io.py`` ではなく ``output.py`` としているのは、
Python標準ライブラリの ``io`` モジュールとの名前衝突を避けるため。
"""


def write_log_header(f):
    f.write("# step  time[fs]  T[K]  KE[Ha]  PE[Ha]  E_tot[Ha]  (E-E0)/|E0|\n")


def write_log_line(f, step, time_fs, T, ke, pe, e_tot, e_drift_rel):
    f.write(
        f"{step:6d}  {time_fs:10.4f}  {T:8.3f}  "
        f"{ke:14.8e}  {pe:14.8e}  {e_tot:14.8e}  {e_drift_rel:+.3e}\n"
    )


def write_xyz_frame(f, pos_bohr, step, bohr_to_angstrom, comment=""):
    """XYZ形式の1フレームを追記する（座標は Å に変換して出力）。"""
    N = pos_bohr.shape[0]
    f.write(f"{N}\n")
    f.write(f"step={step} {comment}\n")
    for i in range(N):
        x = pos_bohr[i, 0] * bohr_to_angstrom
        y = pos_bohr[i, 1] * bohr_to_angstrom
        z = pos_bohr[i, 2] * bohr_to_angstrom
        f.write(f"Ar  {x:12.6f}  {y:12.6f}  {z:12.6f}\n")
