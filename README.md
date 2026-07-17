# chem_3d_builder_rs

2D で描いた化学構造を 3D 化し、分子力学（MM）で構造最適化する統合 GUI アプリケーション。

画面左で 2D 構造式を描画 → 3D 構造へ変換（水素付加＋初期3D座標生成）→ 画面右で力場による最小化を実行する。

## アーキテクチャ

自作の 3 つの crate を接着した上位アプリ:

| crate | 役割 |
|---|---|
| [`chembuider-rs`](https://github.com/IndigoCarmine/chembuider-rs) | 2D 構造エディタ（egui ウィジェット） |
| [`moleucle_3dview_rs`](https://github.com/IndigoCarmine/moleucle_3dview_rs) | 3D 分子ビューア（egui + wgpu） |
| [`openbabel_rs`](https://github.com/IndigoCarmine/openbabel_rs) | 水素付加・3D 生成・MM 力場・最小化・全エクスポート形式 |

エディタで描いた構造は `bridge.rs` で `openbabel::Molecule` へ変換され、以降はそれが唯一の正となる（最小化もエネルギーも書き出しも同じ分子を読む）。

主なソース:

- `src/main.rs` — アプリ本体・GUI レイアウト・エクスポート
- `src/bridge.rs` — 2D → 3D 変換（水素付加・3D 生成）とビューアへの受け渡し
- `src/mm_session.rs` — 最小化の状態管理（ワーカースレッド + 軌跡再生）
- `src/forcefield_kind.rs` — 力場種別の定義

力場は UFF / MMFF94 / MMFF94s / GAFF / Ghemical。OpenBabel 3.2.1 は `forcefieldmm2.cpp` をビルドから除外しているため MM2 は選べない。

## ビルド・実行

GUI は wgpu レンダラを使用（`eframe` の `wgpu` フィーチャ）。egui のバージョンは各サブ crate と合わせて 0.35.x。

`openbabel_rs` は OpenBabel 3.2.1 を**ソースからビルドする**（vcpkg やシステムライブラリは使わない）ため、**CMake と C++ コンパイラが必須**。Windows では MSVC ツールチェインを使うこと。**初回ビルドは 10〜20 分**かかる（2 回目以降は増分）。

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
