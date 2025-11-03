# 指定されたディレクトリ内の全ソースコードを，単一のテキストファイル（.txtまたは.md）に集約するツール
# ソースコードをファイルにまとめる際のファイル名の並び順を，より直感的な「自然順ソート（Natural Sort Order）」に対応させたもの

import os
import sys
import re

try:
    import chardet
except ImportError:
    print("エラー: 'chardet' ライブラリが見つかりません。", file=sys.stderr)
    print(
        "お手数ですが、以下のコマンドでライブラリをインストールしてください:",
        file=sys.stderr,
    )
    print("\npip install chardet\n", file=sys.stderr)
    sys.exit(1)

# --- ここから設定 ---

# [必須] 読み込みたいソースコードが格納されている親フォルダのパス
INPUT_DIR = "/root/wsl/fortran/c_grid/hujiu_yamaken-solver/Code/grid"

# [必須] 生成したマークダウンファイルを出力するフォルダのパス
# このフォルダが存在しない場合は、スクリプトが自動的に作成します。
OUTPUT_DIR = "/root/wsl/fortran/c_grid/hujiu_yamaken-solver/Code"

# [必須] 出力するマークダウンファイル名
OUTPUT_FILENAME = "grid_summarized.md"

# 処理対象とするファイルの拡張子リスト
TARGET_EXTENSIONS = [
    ".f90",  # Fortran
    ".f90inc",  # Fortran include file
    ".c",
    ".h",  # C
    ".cpp",
    ".hpp",
    ".cc",
    ".hh",  # C++
]

# --- ▼▼▼ 変更点: 文字化け対策 ▼▼▼ ---
# ファイルの文字コードを自動判別する際に、優先的に試すエンコーディングのリスト。
# 日本語のコメントが含まれる古いソースコードは 'shift_jis' や 'euc_jp' の可能性があります。
# このリストに挙げられたエンコーディングで読み込みを試み、失敗した場合に chardet による自動判別を行います。
ENCODING_PRIORITY = [
    "utf-8",
    "shift_jis",
    "euc_jp",
    "cp932",  # Shift_JISの亜種
]
# --- ▲▲▲ 変更点 ▲▲▲ ---


# --- 設定はここまで ---


def natural_sort_key(s):
    """
    自然順ソート（例: 'item1', 'item2', 'item10'）のためのキーを生成する。
    """
    return [
        int(text) if text.isdigit() else text.lower()
        for text in re.split(r"([0-9]+)", s)
    ]


def get_language_from_extension(filename):
    """ファイル名から言語名を判定し、Markdownの言語指定子を返す"""
    ext = os.path.splitext(filename)[1].lower()
    if ext in [".f90", ".f90inc"]:
        return "fortran"
    if ext in [".c", ".h"]:
        return "c"
    if ext in [".cpp", ".hpp", ".cc", ".hh"]:
        return "cpp"
    return ""


def append_directory_to_markdown(dir_path, md_file, base_input_dir):
    """
    一つのディレクトリ内のソースファイルを処理し、指定のファイルオブジェクトに追記する
    """
    source_files = []
    # ディレクトリ内の対象ファイルを自然順ソートしてリストアップ
    for filename in sorted(os.listdir(dir_path), key=natural_sort_key):
        ext = os.path.splitext(filename)[1].lower()
        if ext in TARGET_EXTENSIONS:
            file_path = os.path.join(dir_path, filename)
            if os.path.isfile(file_path):
                source_files.append(file_path)

    # 対象ファイルがなければ何もしない
    if not source_files:
        return

    # フォルダ名を見出しとして書き込む
    relative_path = os.path.relpath(dir_path, base_input_dir)
    # ルートディレクトリの場合とサブディレクトリの場合で見出し名を調整
    if relative_path == ".":
        dir_header = os.path.basename(os.path.normpath(base_input_dir))
    else:
        # Windowsのパス区切り文字'\'を'/'に統一
        dir_header = relative_path.replace(os.path.sep, "/")

    md_file.write(f'# folder "{dir_header}"\n\n')
    print(f"\n--- ディレクトリを処理中: {dir_header} ---")

    for file_path in source_files:
        filename = os.path.basename(file_path)
        language = get_language_from_extension(filename)
        print(f"  -> ファイルを追加: {filename}")

        # ファイル内容の読み込み
        content = ""
        try:
            with open(file_path, "rb") as f_binary:
                raw_data = f_binary.read()

            # --- ▼▼▼ 変更点: 文字化け対策の強化 ▼▼▼ ---
            decoded = False
            # 優先エンコーディングリストから順番にデコードを試す
            for encoding in ENCODING_PRIORITY:
                try:
                    content = raw_data.decode(encoding)
                    decoded = True
                    print(f"     - '{filename}' を [{encoding}] で読み込み成功")
                    break
                except UnicodeDecodeError:
                    continue  # デコード失敗、次のエンコーディングへ

            # 優先リストでデコードできなかった場合、chardet を使用
            if not decoded:
                print(
                    f"     - 優先エンコーディングで失敗。chardetで '{filename}' を自動判別します。"
                )
                result = chardet.detect(raw_data)
                encoding = result.get("encoding", "utf-8") or "utf-8"
                confidence = result.get("confidence", "N/A")
                print(f"     - chardet の判別結果: {encoding} (確信度: {confidence})")
                content = raw_data.decode(encoding, errors="replace")
            # --- ▲▲▲ 変更点 ▲▲▲ ---

        except Exception as e:
            content = (
                f"エラー: ファイル '{filename}' の読み込み中に問題が発生しました: {e}"
            )

        # マークダウンへの書き込み (見出しレベルを ## に変更)
        md_file.write(f'## file "{filename}"\n\n')
        md_file.write(f"```{language}\n")
        md_file.write(content.strip())
        md_file.write(f"\n```\n\n")


def main():
    """
    メイン処理
    """
    if not os.path.isdir(INPUT_DIR):
        print(f"エラー: 入力ディレクトリが見つかりません: '{INPUT_DIR}'")
        return

    if not os.path.exists(OUTPUT_DIR):
        print(f"出力ディレクトリ '{OUTPUT_DIR}' を作成します。")
        os.makedirs(OUTPUT_DIR)

    output_md_path = os.path.join(OUTPUT_DIR, OUTPUT_FILENAME)

    print("処理を開始します...")
    print(f"入力ディレクトリ: {os.path.abspath(INPUT_DIR)}")
    print(f"出力ファイル: {os.path.abspath(output_md_path)}")

    try:
        with open(output_md_path, "w", encoding="utf-8") as md_file:
            # os.walkでINPUT_DIR以下の全サブディレクトリを探索
            # 処理順序を固定するため、ディレクトリパスでソート
            sorted_walk = sorted(
                os.walk(INPUT_DIR), key=lambda x: natural_sort_key(x[0])
            )

            for dir_path, _, _ in sorted_walk:
                append_directory_to_markdown(dir_path, md_file, INPUT_DIR)

        print("\nすべての処理が完了しました。")

    except IOError as e:
        print(
            f"エラー: ファイル '{output_md_path}' の書き込み中に問題が発生しました: {e}",
            file=sys.stderr,
        )


if __name__ == "__main__":
    main()
