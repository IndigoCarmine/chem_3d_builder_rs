# chem_3d_builder_rs

2D で描いた化学構造を 3D 化し、分子力学（MM）で構造最適化する統合 GUI アプリケーション。

画面左で 2D 構造式を描画 → 3D 構造へ変換（水素付加＋初期3D座標生成）→ 画面右で力場による最小化を実行する。

## アーキテクチャ

自作の 3 つの crate を接着した上位アプリ:

| crate | 役割 |
|---|---|
| [`chembuider-rs`](https://github.com/IndigoCarmine/chembuider-rs) | 2D 構造エディタ（egui ウィジェット） |
| [`moleucle_3dview_rs`](https://github.com/IndigoCarmine/moleucle_3dview_rs) | 3D 分子ビューア（egui + wgpu） |
| [`openbabel_rs`](https://github.com/IndigoCarmine/openbabel_rs) | 水素付加・3D 生成・MM 力場・最小化・二面角・全入出力形式 |

エディタで描いた構造は `bridge.rs` で `openbabel::Molecule` へ変換され、以降はそれが唯一の正となる（最小化もエネルギーも書き出しも同じ分子を読む）。

主なソース:

- `src/main.rs` — アプリ本体・GUI レイアウト・ファイル入出力
- `src/bridge.rs` — 2D → 3D 変換（水素付加・3D 生成）とビューアへの受け渡し
- `src/mm_session.rs` — 最小化の状態管理（ワーカースレッド + 軌跡再生）
- `src/forcefield_kind.rs` — 力場種別の定義
- `src/viewport3d.rs` — 3D ビューポート（後述の理由で自前実装）
- `src/edit3d.rs` — 3D 編集の選択状態と右パネル
- `src/dihedral.rs` — 二面角の対話回転（順運動学）
- `src/ik.rs` — 2 原子間距離から二面角を解く（逆運動学）
- `src/geom3d.rs` — 回転・グラフ探索・平面判定・原子の重なり検出
- `src/test_support.rs` — OpenBabel を叩くテストの直列化（テストビルドのみ）

3D ビューポートは `moleucle_3dview_rs` の `InteractiveMoleculeViewport` を使わず、
同 crate の公開部品（`MoleculeViewer` / `OrbitalCamera` / `OffscreenRenderer` /
`RenderFrameState`）から `src/viewport3d.rs` で組み立てている。上流のウィジェットは
原子の左クリックしか外に出さず（`pick()` が返す `BondClicked` は内部で捨てられる）、
ホイールはカメラ、右ボタンはパンに直結しているため、下記の 3D 編集操作が実装できない。
描画部分は上流の `show()` をそのまま写しており、入力処理だけが異なる（`DIFFERS` コメント）。

## 3D 生成の信頼性

OpenBabel の 3D ビルダーは確率的で、同じ構造でも 8 回に 1 回ほど壊れた配置を返す。
`bridge.rs` は最大 5 回試し、結果を 3 つに分類する:

| 判定 | 条件 | 扱い |
|---|---|---|
| Clean | 問題なし | そのまま採用 |
| Cramped | 非結合の原子対が 1.2 Å 未満 | **警告付きで採用**（MM で解消できる）|
| Unusable | NaN 座標 / 座標が 1e4 Å 超へ飛んだ / 実質同一点 / 座標が読めない | 破棄して再試行 |

さらに、採用する前に**選択中の力場でエネルギーを計算し、有限かつ 1e5 以下**であることを確認する。
距離だけでは足りないため:

- ビルダーは NaN ではなく**有限だが天文学的な座標**（実測 1e243 Å）を書くことがある。
  1 原子が飛ぶと残りが 0.15 Å に潰れるので、距離だけ見ると「ちょっと近い」に見えてしまう。
- 距離が正常でもエネルギーが 1e14 になる配置が実測で 1.4% あった。
  ユーザーが見るのはエネルギーなので、それ自体を合否の基準にする。

400 回の生成で拒否 0 件・異常エネルギー 0 件（対策前は最悪 3.2e14）。
なお `Badge` は非有限の値では作られないので、`E NaN` が表示されることはない。

## 3D 編集

生成した立体構造は 3D ビュー上で直接いじれる。右パネル「3D 編集」が選択状態と操作を持つ。

| 操作 | 効果 |
|---|---|
| 原子を右クリック | アンカー原子（動かない側の目印）にする |
| 結合を左クリック | 二面角の回転軸にする |
| ホイール | アンカーが属さない側だけを回す（1 ノッチ 5°、Shift で 1°、Ctrl で 15°）|
| 原子を左クリック | 逆運動学の対象原子にする（2 つまで、3 つ目で古い方を捨てる）|
| 何もない所を左/右クリック | 結合と対象原子／アンカーの選択を解除 |

アンカーと結合が揃っているときだけホイールがアプリ側に渡り、それ以外はこれまで通り
カメラのズームになる。環内・多重結合・末端の結合は `Bond::is_rotor()` が false なので
拒否し、理由をパネルに出す。

**逆運動学**は 2 原子を選び、目標（距離を指定 / できるだけ遠く / できるだけ近く）を与えると、
その 2 原子の最短経路上にある回転可能結合の二面角を解く。1 自由度あたり距離の 2 乗が
`d²(φ) = C − 2R·cos(φ − α)` の閉形式になるので、CCD（各結合を順に閉形式で最適化）で
数回の反復で収束する。指定距離に届かない場合は到達可能な範囲を添えて報告し、
原子の重なりが生じた場合は件数と最悪の対を警告する（MM は自動では掛けない）。

3D 編集は直前 20 手ぶんの座標を保持しており、「↶ 3D 編集を取り消す」で戻せる
（MM 最小化を回すとこの履歴は破棄される）。2D エディタの undo は 2D 構造だけを見るので、
両者は独立している。

## ファイル入出力

| 操作 | 形式 |
|---|---|
| ⭱ 読み込み (Ctrl+O) | MDL Molfile / SDF、Mol2、PDB、XYZ、SMILES、InChI |
| ⭳ エクスポート (Ctrl+S) | Mol2、PDB、SDF、XYZ、SMILES |
| 🖼 画像 | 3D ビュー (PNG)、2D 構造式 (SVG) |

形式は拡張子から OpenBabel が推定する。読み書きの対応表は
`READ_FORMATS` / `FORMATS`（`src/main.rs`）にあり、`reads_every_offered_format` /
`writes_every_offered_format` が「ダイアログに出す拡張子はすべて実際に読み書きできる」ことを
テストで固定している。読みと書きは別のプラグインなので、両方を別々に確認している。
CML と CIF は OpenBabel が libxml2 に対してビルドされていないため出していない。

読み込みには 2 つ制約がある:

- **左の 2D キャンバスには反映されない。** `bridge.rs` は 2D → OpenBabel の片道変換しか
  持たないため、読み込んだ構造を 2D で描き直す手段がない。読み込み後に「→ 3D 生成」を
  押すと、読み込んだ構造は 2D の内容で上書きされる。
- **平面的な構造は立体構造を生成し直す。** Mol2 や SDF は座標が平面でも次元を 3 と申告する
  ため、`has_3d()` だけでは 2D の下書きを見分けられない（そのまま読むと原子が重なり、
  エネルギーが数億 kJ/mol になる）。座標そのものが平面かどうかを見て、平面なら水素付加と
  3D 生成をやり直し、その旨をバッジに出す。

力場は UFF / MMFF94 / MMFF94s / GAFF / Ghemical。OpenBabel 3.2.1 は `forcefieldmm2.cpp` をビルドから除外しているため MM2 は選べない。

## ビルド・実行

GUI は wgpu レンダラを使用（`eframe` の `wgpu` フィーチャ）。egui のバージョンは各サブ crate と合わせて 0.35.x。

`openbabel_rs` は OpenBabel 3.2.1 を**ソースからビルドする**（vcpkg やシステムライブラリは使わない）ため、**CMake と C++ コンパイラが必須**。Windows では MSVC ツールチェインを使うこと。**初回ビルドは 10〜20 分**かかる（2 回目以降は増分）。
CI（と Windows インストーラーのビルド）では、`openbabel_rs` のリリースに添付されたビルド済み OpenBabel を `OPENBABEL_SYS_PREBUILT_DIR` で使ってこのビルドを省略している。そのため `Cargo.toml` の `openbabel` は**リリース済みのタグ**で固定すること（詳細は openbabel_rs の [Prebuilt OpenBabel for CI](https://github.com/IndigoCarmine/openbabel_rs/blob/master/docs/src/building.md#prebuilt-openbabel)）。

```sh
cargo run --release
```

## リリース / インストーラー

バージョンタグ（`v*`、例 `v0.4.0`）を push すると、GitHub Actions（[`.github/workflows/release.yml`](.github/workflows/release.yml)）が 3 OS 分のネイティブインストーラーをビルドし、同名の GitHub Release に添付する。

| OS | 生成物 | ツール |
|---|---|---|
| Windows | `chem_3d_builder-<ver>-setup.exe` | NSIS（[`packaging/windows/installer.nsi`](packaging/windows/installer.nsi)） |
| macOS | `chem_3d_builder-<ver>.dmg`（ユニバーサルバイナリ） | `.app` バンドル + `hdiutil` |
| Linux | `chem_3d_builder-<ver>-x86_64.AppImage` | `linuxdeploy` + gtk プラグイン |

リリース手順:

```sh
# Cargo.toml の version を上げてコミットしてから
git tag v0.4.0
git push origin v0.4.0
```

アイコンは [`assets/icon.svg`](assets/icon.svg) を各形式（PNG / ICNS）へ CI 内で変換して使用する（Windows は既定アイコン）。差し替えたい場合はこの SVG を編集する。

> インストーラーはコード署名していないため、初回起動時に SmartScreen（Windows）や Gatekeeper（macOS: 右クリック →「開く」）の警告が出る。署名するには各社の証明書とワークフローへの署名ステップ追加が必要。

### OpenBabel ランタイムの同梱

インストーラーは実行ファイルの隣に OpenBabel のランタイムを並べる。このレイアウト自体が動作条件で、飾りではない:

```text
chem_3d_builder_rs(.exe)
openbabel-3.dll / libopenbabel.*   共有ライブラリ
*.obf                              フォーマット・力場プラグイン
data/                              力場の .par 等
```

`data/` が無いとアプリは**普通に起動したうえで**力場と 3D 生成だけが黙って死ぬ（`Cannot open UFF.prm` / `OBBuilder::LoadFragments`）。ビルドツリーには焼き込みパスがあるためローカルでは絶対に再現せず、配布版でのみ起きる。`.obf` は `dlopen` されるので linuxdeploy のような依存追跡ツールも見つけられない。ワークフローの各パッケージングステップは、この 3 つが揃わなければ Release を失敗させる。

## ライセンス

GPL-2.0-only。OpenBabel（GPL-2.0-only）をリンクするため、結合著作物がこのライセンスを継承する。
