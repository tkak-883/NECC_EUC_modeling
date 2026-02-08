# NECC / EUC 再現モデル (Fortran + OpenMP)

北太平洋の赤道帯を **理想化した水槽**として扱い、有限差分法（リープフロッグ + Asselin フィルター）で
**北赤道反流 (NECC)** と **赤道潜流 (EUC)** の再現を試みる Fortran コードです。
計算コスト削減のため **OpenMP 並列化**も入っています。

- NECC: 表層で北緯 ~15° 付近に東向き（~0.2 m/s 程度）
- EUC : 深さ ~100–200 m 付近に東向きジェット（~0.3 m/s 程度）
  （※パラメータによって変化）

---

## 特徴

- 3次元プリミティブ方程式型の簡易モデル (ブシネスク + 静水圧近似)
- 予報変数: `u, v, eta (水位), T, S`
- 診断変数: `w, rho, p`
- 風応力強制: 貿易風 + 偏西風 (Gauss分布)
- 鉛直対流調整
- OpenMPによる並列化(`collapse`, `schedule(static)`, `workshare`, `single`)

---

## 実行環境

- `gfortran` (OpenMP対応)
- macOS / Linux を想定（Windowsでもgfortranが使えれば可）
- OpenMP: `-fopenmp`

---

## ビルド

```bash
gfortran -O3 -fopenmp report_parallel_necc.f90 -o report_parallel_necc
```

※ファイル名は適宜置き換えてください。

---

## 実行

スレッド数を環境変数で指定して実行します。

```bash
export OMP_NUM_THREADS=4
./report_parallel_necc
```

スレッド数の例：

```bash
export OMP_NUM_THREADS=1   # 並列化オーバヘッドあり
export OMP_NUM_THREADS=4   # この環境では最速になりやすい
export OMP_NUM_THREADS=8   # M2のEコア混在で遅くなる場合あり
```

---

## モデル / パラメータ (default)

- グリッド数: `im=100, jm=40, km=12`
- グリッド幅: `dx=100 km, dy=100 km, dz=50 m` → depth = 600 m
- 時間幅: `dt=300 s`
- 計算期間: `5 years`
- Output interval: `30 days`
- コリオリ力: **β-plane のみ (f = βy)**
- 粘性 / 拡散係数:
  - 水平粘性係数 `ah=1e5 m^2/s`
  - 鉛直粘性係数 `av=1e-3 m^2/s`
  - 鉛直拡散係数 `kv=1e-4 m^2/s`
  - 海面高度の水平拡散係数 `kh=0`（今回は使用せず）

風応力 (example):

```fortran
tau_x(y) = -tau_tr * exp(-((lat-lat_tr)/wid_tr)^2) + tau_w * exp(-((lat-lat_w)/wid_w)^2)
```

---

## Output files

`expid`（デフォルト: `ex4-case1`）を接頭辞として以下を出力します。

- `ex4-case1.nXXXXXX` : 解析用スナップショット（unformatted, stream）
  - `time_sec, u, v, eta, w, p, T, S, rho`

- `ex4-case1.cnt` : 継続計算用チェックポイント（unformatted, stream）

> **注意**: `unformatted` + `access='stream'` のため、読み出しは Fortran / Python 側で dtype を合わせてください。

---

## 可視化

手元の可視化スクリプト例（別途用意）：

- `ex4.py` : 水平面（流速ベクトル、u など）
- `ex4_vertical_section.py` : 南北-鉛直断面（u, 温度, 塩分など）

（例：経度 `i=50`、水平断面 `k=9,12` など）

---

## スレッド数について

このコードは 3重ループが支配的なので OpenMP が効きますが、
Apple Silicon (M2) は **Pコア + Eコア混在**のため、スレッド数を増やしすぎると逆に遅くなることがあります。

実測例（同条件）：

- original: 6m17s
- output interval = 5 years: 6m12s
- division → multiplication (precompute reciprocals): 5m31s
- OpenMP 4 threads: **4m46s**
- OpenMP 8 threads: 7m11s
- OpenMP 1 thread: 6m11s

また、並列化後も同一条件で **結果が一致**することを確認済みです。

---

## Acknowledgements

授業で配布されたサンプルコード・講義資料を参考にしつつ、
風応力分布（NECC/EUC向け）と OpenMP 並列化を追加して実験したものです。
