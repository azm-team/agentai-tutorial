#import "@preview/cetz:0.4.2"
#import "@preview/cetz-plot:0.1.3": plot

// ----- ドキュメント設定 -----
#set document(title: "分子動力学（MD）入門", author: "agentai-tutorial")
#set page(
  paper: "a4",
  margin: (x: 2cm, top: 1.8cm, bottom: 2cm),
  numbering: "1 / 1",
  number-align: center,
)
#set text(
  font: ("Hiragino Mincho ProN", "Noto Serif CJK JP"),
  size: 10pt,
  lang: "ja",
)
#set par(justify: true, leading: 0.7em, first-line-indent: 1em)
#set heading(numbering: "1.1")
#show heading.where(level: 1): it => block(above: 1.2em, below: 0.6em)[
  #set text(size: 13pt, weight: "bold")
  #it
]
#show heading.where(level: 2): it => block(above: 0.9em, below: 0.4em)[
  #set text(size: 11pt, weight: "bold")
  #it
]
#show figure.caption: set text(size: 9pt)
#show raw: set text(font: "Menlo", size: 9pt)
#show link: set text(fill: rgb("0050b0"))

// ----- 表紙 -----
#align(center)[
  #v(0.2em)
  #text(size: 18pt, weight: "bold")[分子動力学（MD）入門]
  #v(0.3em)
  #text(size: 12pt)[―― 理論的基礎と便利なツール ――]
  #v(0.3em)
  #text(size: 9pt)[agentai-tutorial / 講義ノート]
]
#v(0.5em)

= はじめに

分子動力学（Molecular Dynamics, MD）は，原子核を古典粒子とみなしてその運動方程式を数値的に解くシミュレーション手法である。原子間に働く力 $bold(F) = - nabla V$ から各原子の加速度を求め，時間刻み $Delta t$ ずつ位置と速度を更新する。これにより，

- 原子配置の時間発展（液体・固体の構造，自己拡散，相転移）
- 熱力学量（温度，圧力，弾性率，比熱）
- 自由エネルギー（FEP, MetaD, Umbrella Sampling 等）

などが計算できる。本稿では (i) MDを理解するうえで必要な理論を最小限まとめ，続いて (ii) 実際に使われる主要ツールを概観する。

= MDの理論的基礎

== 古典力学の復習

$N$ 個の原子からなる系のハミルトニアンは

$ H(bold(r), bold(p)) = sum_(i=1)^N frac(|bold(p)_i|^2, 2 m_i) + V(bold(r)_1, dots, bold(r)_N) $

で与えられる。各原子の運動方程式は Newton の第二法則

$ m_i frac(d^2 bold(r)_i, d t^2) = - nabla_i V $

であり，これを数値積分して軌跡 $bold(r)_i (t)$ を得る。MDの主目的はこの軌跡から熱力学的平均量を取り出すこと（ergodic仮説のもと，時間平均 ≈ アンサンブル平均）にある。

== ポテンシャルエネルギー $V$ の表現

最も単純な例が Lennard-Jones (LJ) ポテンシャル：

$ V_(L J) (r) = 4 epsilon [(sigma / r)^12 - (sigma / r)^6] $

第1項は短距離反発（Pauli斥力），第2項は分散引力（London力）を表す。$sigma$ は $V = 0$ となる距離，$epsilon$ はポテンシャル井戸の深さで，最小は $r_min = 2^(1/6) sigma approx 1.122 sigma$ にある。希ガス（Ar, Kr等）の良い近似となる。

#figure(
  cetz.canvas({
    import cetz.draw: *
    plot.plot(
      size: (10, 4.5),
      x-label: $r \/ sigma$,
      y-label: $V \/ epsilon$,
      x-min: 0.9, x-max: 3.0,
      y-min: -1.4, y-max: 2.0,
      x-tick-step: 0.5,
      y-tick-step: 1.0,
      axis-style: "school-book",
      {
        plot.add(
          domain: (0.95, 3.0),
          samples: 200,
          x => 4 * (calc.pow(x, -12) - calc.pow(x, -6)),
          style: (stroke: 1.5pt + blue),
          label: $V_(L J)$,
        )
        plot.add(
          ((calc.pow(2, 1.0/6.0), -1.0),),
          mark: "o", mark-size: 0.2, mark-style: (fill: red, stroke: red),
        )
      },
    )
  }),
  caption: [Lennard-Jonesポテンシャル。最小値は $r = 2^(1/6) sigma$（赤丸）で $V = -epsilon$。$r < sigma$ で急峻に発散，$r → infinity$ で 0 に漸近する。],
)

実用上は，生体・有機系では結合伸縮（bond）・角度（angle）・二面角（torsion）・非結合相互作用（vdW + Coulomb）を組み合わせた力場が用いられる：

$ V = V_"bond" + V_"angle" + V_"torsion" + V_(L J) + V_"Coulomb" $

これらの関数形と係数を集めたものが AMBER, CHARMM, OPLS などの「力場（force field）」である。

== 時間積分：Velocity Verlet法

MDで最も広く使われる積分法が *Velocity Verlet* 法：

$ bold(v)(t + Delta t \/ 2) &= bold(v)(t) + frac(Delta t, 2 m) bold(F)(t) quad &&"(half-kick)" \
  bold(r)(t + Delta t) &= bold(r)(t) + Delta t bold(v)(t + Delta t \/ 2) quad &&"(drift)" \
  bold(F)(t + Delta t) &= - nabla V(bold(r)(t + Delta t)) quad &&"(force update)" \
  bold(v)(t + Delta t) &= bold(v)(t + Delta t \/ 2) + frac(Delta t, 2 m) bold(F)(t + Delta t) quad &&"(half-kick)" $

時間反転対称・シンプレクティック性を持ち，NVEで長時間にわたりエネルギー保存が良い（局所誤差 $O(Delta t^3)$，全エネルギー誤差 $O(Delta t^2)$ で振動）。

#figure(
  cetz.canvas({
    import cetz.draw: *
    let box(pos, label) = {
      rect(
        (pos.at(0) - 1.4, pos.at(1) - 0.45),
        (pos.at(0) + 1.4, pos.at(1) + 0.45),
        radius: 0.1, stroke: 0.8pt, fill: rgb("eef"),
      )
      content(pos, text(size: 8pt, label))
    }
    box((0, 0), [half-kick \ v ← v + (dt/2)F/m])
    box((3.5, 0), [drift \ r ← r + dt·v])
    box((7, 0), [force update \ F ← −∇V(r)])
    box((10.5, 0), [half-kick \ v ← v + (dt/2)F/m])
    line((1.4, 0), (2.1, 0), mark: (end: ">"))
    line((4.9, 0), (5.6, 0), mark: (end: ">"))
    line((8.4, 0), (9.1, 0), mark: (end: ">"))
    bezier(
      (11.9, 0.45), (-1.4, 0.45),
      (12.5, 1.8), (-2.0, 1.8),
      stroke: 0.6pt + gray, mark: (end: ">"),
    )
    content((5.25, 1.7), text(size: 8pt, fill: gray)[次のステップへ])
  }),
  caption: [Velocity Verlet 1ステップ（kick-drift-force-kick）。1関数にまとめると順序ミスを防げる。],
)

== 周期境界条件と最小イメージ規約

有限サイズ系の表面効果を避けるため，立方体ボックスを無限に並べた周期境界条件（PBC）を採る。粒子間距離は最も近い周期像との距離（minimum image）で測る：

$ Delta r_alpha "→" Delta r_alpha - L dot "round"(Delta r_alpha / L), quad alpha in {x, y, z} $

カットオフ距離 $r_c$ は $L \/ 2$ より小さく取らねばならない。さらに，$r_c$ でのポテンシャル不連続を避けるため $U(r_c)$ を引き算する *energy shift* を施す（NVEエネルギー保存に必須）。

#figure(
  cetz.canvas({
    import cetz.draw: *
    let atoms = ((0.5, 0.4), (1.5, 0.7), (1.3, 1.6), (0.4, 1.4))
    for ix in (-1, 0, 1) {
      for iy in (-1, 0, 1) {
        let is_center = (ix == 0 and iy == 0)
        let s = if is_center { 1pt + black } else { 0.4pt + gray }
        let f = if is_center { rgb("eef") } else { none }
        rect((ix * 2, iy * 2), (ix * 2 + 2, iy * 2 + 2), stroke: s, fill: f)
        for a in atoms {
          let r = if is_center { 0.15 } else { 0.10 }
          let c = if is_center { blue } else { rgb("88a") }
          circle((a.at(0) + ix * 2, a.at(1) + iy * 2), radius: r, fill: c, stroke: none)
        }
      }
    }
    // 最小イメージベクトル例：中央セルの atom0 から左セルの atom2 像へ
    let from = atoms.at(0)
    let to = (atoms.at(2).at(0) - 2, atoms.at(2).at(1))
    line(from, to, stroke: 1pt + red, mark: (end: ">"))
    content((-0.4, 2.6), text(size: 8pt, fill: red)[最小イメージ])
    content((1, -1.5), text(size: 8pt)[中央：シミュレーションボックス（実体），周囲：周期像])
  }),
  caption: [2次元 PBC の模式図。各ペアは最も近い周期像との距離で相互作用を計算する。],
)

== アンサンブル

#figure(
  table(
    columns: (auto, auto, 1fr),
    align: (center, center, left),
    stroke: 0.5pt,
    [*略号*], [*保存量*], [*必要な機構*],
    [NVE], [$N, V, E$], [Velocity Verlet のみで実現（最も基本）],
    [NVT], [$N, V, T$], [サーモスタット：Berendsen / Nosé–Hoover / Langevin / VRescale],
    [NPT], [$N, P, T$], [サーモ＋バロスタット：Berendsen / Parrinello–Rahman / MTK],
  ),
  caption: [代表的なアンサンブルと実現機構。],
)

実験条件に近いのは NPT。NVE はエネルギー保存の検証や輸送物性計算（自己拡散係数等）で使う。

== 初期配置と初期速度

固体・液体の出発点として FCC 格子（面心立方）が便利。1 ユニットセルに 4 原子の basis を持つ：

$ (0,0,0), space (1/2, 1/2, 0), space (1/2, 0, 1/2), space (0, 1/2, 1/2) $

各原子に Maxwell–Boltzmann 分布から速度を与え（$sigma_v = sqrt(k_B T \/ m)$），重心並進を除去後，目標温度 $T_0$ に厳密リスケールする（自由度は $3N - 3$，COM を除いた数）。

#figure(
  cetz.canvas({
    import cetz.draw: *
    rect((0, 0), (3.5, 3.5), stroke: 0.8pt)
    // 角（z=0, z=a）
    for (x, y) in ((0, 0), (3.5, 0), (0, 3.5), (3.5, 3.5)) {
      circle((x, y), radius: 0.28, fill: blue)
    }
    // 面心 (1/2,1/2,0)
    circle((1.75, 1.75), radius: 0.28, fill: red)
    // 面心 (1/2,0,1/2) と (1/2,1,1/2) → 投影で重なる
    circle((1.75, 0), radius: 0.28, fill: rgb("00a000"))
    circle((1.75, 3.5), radius: 0.28, fill: rgb("00a000"))
    // 面心 (0,1/2,1/2) と (1,1/2,1/2)
    circle((0, 1.75), radius: 0.28, fill: orange)
    circle((3.5, 1.75), radius: 0.28, fill: orange)
    // 凡例
    let lg(y, c, label) = {
      circle((5.0, y), radius: 0.18, fill: c)
      content((5.4, y), align(left, text(size: 9pt, label)))
    }
    lg(3.0, blue, [角 $(0,0,0)$ 系])
    lg(2.3, red, [面心 $(1/2,1/2,0)$])
    lg(1.6, rgb("00a000"), [面心 $(1/2,0,1/2)$])
    lg(0.9, orange, [面心 $(0,1/2,1/2)$])
    content((1.75, -0.6), text(size: 8pt)[FCC unit cell ($x y$ 投影)])
  }),
  caption: [FCC ユニットセルの $x y$ 平面投影。角と各面の中心に原子があり，1セルあたり実質 4 原子。],
)

= 便利なツール

== MDエンジン（汎用シミュレーションコード）

#figure(
  table(
    columns: (auto, 1fr, auto, 2fr),
    align: (left, left, center, left),
    stroke: 0.5pt,
    [*エンジン*], [*主用途*], [*GPU*], [*特徴*],
    [LAMMPS], [一般／材料系（金属・酸化物・ポリマー）], [○], [プラグイン豊富，スクリプト言語で柔軟。MLIP統合多数。],
    [GROMACS], [生体高分子・脂質膜], [○], [生体MDで最速級。AMBER/CHARMM/OPLS対応。],
    [NAMD], [大規模生体系], [○], [並列スケーラビリティ重視。VMDと密連携。],
    [OpenMM], [研究・教育・MLIP], [○], [Python API。MLIPやカスタム力との統合が容易。],
  ),
  caption: [主要MDエンジンの比較。],
)

== Pythonライブラリ

- *ASE* (Atomic Simulation Environment) ── 構造構築・各種計算機（VASP, QE, LAMMPS, MLIP）の統一インターフェイス。MD・構造最適化・NEB も実装。
- *MDAnalysis* ── 軌跡解析（RDF, MSD, RMSD等）の事実上の標準。多数のフォーマットを読み込み可能。
- *MDTraj* ── 軌跡解析・フォーマット変換。一部関数で GPU 加速あり。

== 可視化ツール

- *VMD* ── 生体系の定番。Tcl/Pythonで自動化可能。表現が豊富。
- *OVITO* ── 結晶欠陥・粒界・転位解析に強い。Python モジュールも提供。
- *PyMOL* ── タンパク質可視化に強い。論文用Figure作成によく使われる。

== 力場・ポテンシャル

*古典力場*：

- *AMBER, CHARMM* ── タンパク質・核酸・脂質
- *OPLS-AA / GAFF* ── 有機分子・溶媒
- *Lennard-Jones, Morse* ── 希ガス，単純流体（教育用にも）
- *EAM, MEAM* ── 金属系

特徴：軽量（$mu s$ 規模のシミュレーション可），だが反応や電子状態変化は扱えない。

*機械学習ポテンシャル（MLP / MLIP）*：

従来の NNP（Behler–Parrinello, GAP, ANI など）は系ごとに DFT 訓練データが必要だった。近年は大規模事前学習による *Universal MLIP*（基盤MLIP）が登場し，任意元素系に転用可能となっている：

#figure(
  table(
    columns: (auto, auto, 1fr),
    align: (left, left, left),
    stroke: 0.5pt,
    [*モデル*], [*対象元素*], [*特徴*],
    [MACE-MP-0], [周期表ほぼ全域], [Equivariant message-passing（高次テンソル）。Materials Project + GNoME 等で事前学習。],
    [MatterSim], [元素・温度・圧力広範], [Microsoft Research。条件付き訓練で広いP-T範囲を網羅。],
    [M3GNet], [89 元素], [構造-エネルギー-力-応力を統一的に予測。],
    [CHGNet], [周期表], [磁気モーメントも同時予測。],
    [SevenNet], [周期表], [Equivariant transformer。並列推論最適化。],
  ),
  caption: [代表的な Universal MLIP。],
)

利点：DFTに近い精度（典型誤差 数 meV/atom 程度）で，事前学習済みモデルを *zero-shot* 適用可能。古典力場よりは1〜3桁重いが，DFT MD よりは $10^3 tilde.op 10^5$ 倍速い。新材料探索・触媒・電解液・固体電解質などで急速に普及している。

ASE や OpenMM, LAMMPS から呼び出すラッパが整備されており，導入の敷居は急速に下がっている。

= 入門者向け Tips

- *時間刻み*：最も速い振動周期の 1/10 以下。LJ-Ar なら 1–5 fs。生体系では X–H 結合振動を SHAKE/LINCS で固定して 2 fs。
- *カットオフ距離*：ボックス長の半分以下。LJ は energy shift か smooth truncation で連続化。Coulomb は Ewald / PME を使う。
- *NVE 検証*：開発時はまず NVE で $|Delta E \/ E_0|$ を確認する。これが満たせないなら $Delta t$，カットオフ，力計算のいずれかが怪しい。
- *単位系*：reduced (LJ) / real (Å, fs, kcal/mol) / metal (Å, ps, eV) / atomic (bohr, Hartree) など複数ある。入力ファイル冒頭で必ず明示・確認する。
- *並列化*：最初は 1 ノード・短ステップで動かして物理を確認。GPU/MPI スケールはその後。
- *再現性*：初期速度の乱数 seed と全パラメータを記録する。MDは確率的なので「実行のたびに違う結果」は当然。

= 参考文献

- D. Frenkel and B. Smit, _Understanding Molecular Simulation_, 3rd ed., Academic Press (2023).
- M. P. Allen and D. J. Tildesley, _Computer Simulation of Liquids_, 2nd ed., Oxford (2017).
- LAMMPS Documentation: #link("https://docs.lammps.org/")
- GROMACS Documentation: #link("https://manual.gromacs.org/")
- OpenMM Documentation: #link("https://docs.openmm.org/")
- I. Batatia et al., "MACE: Higher Order Equivariant Message Passing Neural Networks for Fast and Accurate Force Fields," _NeurIPS_ (2022).
- H. Yang et al., "MatterSim: A Deep Learning Atomistic Model Across Elements, Temperatures and Pressures," _arXiv:2405.04967_ (2024).
