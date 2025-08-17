#!/bin/bash

cd ~/wsl/python/sudoku-camera

# 1-9のテンプレートを作成
python shogicam/preprocess/make_templates_stronger.py

# ラベルを手動で設定するために81個の画像を収集
python shogicam/preprocess/harvest_unlabeled.py

# テンプレートのラベルを手動で設定し、それを自動でテンプレートに保存
python shogicam/preprocess/build_templates_from_labeled.py

# テンプレートを元に画像認識
python -m shogicam.preprocess.board_detection