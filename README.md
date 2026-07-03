# chem_3d_builder_rs

2D で描いた化学構造を 3D 化し、分子力学（MM）で構造最適化する統合 GUI アプリケーション。

画面左で 2D 構造式を描画 → 3D 構造へ変換（水素付加＋初期3D座標生成）→ 画面右で力場による最小化を実行する。

## アーキテクチャ

自作の 3 つの crate を接着した上位アプリ:

| crate | 役割 |
|---|---|
| [`chembuider-rs`](https://github.com/IndigoCarmine/chembuider-rs) | 2D 構造エディタ（egui ウィジェット） |
| [`moleucle_3dview_rs`](https://github.com/IndigoCarmine/moleucle_3dview_rs) | 3D 分子ビューア（egui + wgpu） |
| [`zunda_rs`](https://github.com/IndigoCarmine/zunda_rs) | MM 力場・最小化（UFF / Ghemical 等） |

主なソース:

- `src/main.rs` — アプリ本体・GUI レイアウト
- `src/bridge.rs` — 2D → 3D 変換（水素付加・初期3D座標）
- `src/mm_session.rs` — MM 最小化セッションの管理
- `src/forcefield_kind.rs` — 力場種別の定義

## ビルド・実行

GUI は wgpu レンダラを使用（`eframe` の `wgpu` フィーチャ）。egui のバージョンは各サブ crate と合わせて 0.34.x。

```sh
cargo run --release
```
