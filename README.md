# ElasticWaveFDTD-CUDA

**ElasticWaveFDTD-CUDA**は、弾性波の有限差分時間領域法（FDTD）をCUDA C++で実装した高性能なシミュレーションライブラリです。本ライブラリは、研究者やエンジニアが地震波のシミュレーションや材料科学、医療用途など多岐にわたる応用分野で利用できるよう設計されています。

## 目次
- [特徴](#特徴)
- [デモ](#デモ)
- [インストール](#インストール)
- [使用方法](#使用方法)
  - [基本的な使い方](#基本的な使い方)
  - [サンプルプロジェクト](#サンプルプロジェクト)
- [APIドキュメント](#apiドキュメント)
- [パフォーマンス](#パフォーマンス)
- [拡張機能](#拡張機能)
- [貢献方法](#貢献方法)
- [ライセンス](#ライセンス)
- [お問い合わせ](#お問い合わせ)
- [参考文献](#参考文献)

## 特徴

- **高性能なCUDA実装**: NVIDIA GPUを活用し、計算速度を最大化。
- **適応メッシュ細分化（AMR）**: シミュレーション領域の特定部分でメッシュを動的に細分化し、精度と効率を両立。
- **高次精度スキーム**: 2次および4次の有限差分スキームをサポートし、高精度な結果を提供。
- **非線形材料特性のモデル化**: 線形および非線形弾性材料の特性をシミュレーション可能。
- **シンプルで直感的なAPI**: ユーザーフレンドリーなインターフェースにより、容易なシミュレーション設定と実行を実現。
- **マルチGPUサポート**: 複数のGPUを活用した並列計算に対応。
- **オープンソース**: GitHub上で公開され、コミュニティによる拡張と改善が可能。

## デモ

以下のGIFは、ElasticWaveFDTD-CUDAを使用してシミュレーションされた弾性波の伝搬を示しています。

![デモ映像](https://example.com/demo.gif)

## インストール

### 前提条件

- **CUDA対応GPU**: NVIDIA GPUが必要です。
- **CUDA Toolkit**: バージョン11.0以上を推奨。
- **C++コンパイラ**: C++17対応のコンパイラ（例: `gcc`, `clang`）。
- **CMake**: バージョン3.10以上。

### 手順

1. **リポジトリをクローン**

   ```bash
   git clone https://github.com/yourusername/ElasticWaveFDTD-CUDA.git
   cd ElasticWaveFDTD-CUDA