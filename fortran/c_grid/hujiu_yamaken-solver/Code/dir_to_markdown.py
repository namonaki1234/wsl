"""
dir_to_markdown_fixed.py
指定したディレクトリを再帰的に走査して Markdown のツリーを出力する
"""

import os
from pathlib import Path

# ======= 設定セクション =======
# ここを書き換えるだけでOK
ROOT_DIR = Path("/root/wsl/fortran/c_grid/hujiu_yamaken-solver/Code/solver")
OUT_PATH = Path("/root/wsl/fortran/c_grid/hujiu_yamaken-solver/Code/solver_tree.md")
# =================================


def make_md_tree_from_fs(root: Path) -> str:
    """ディレクトリを再帰的に Markdown ツリー化"""
    root = root.resolve()
    lines = []
    lines.append("```text")
    lines.append(f"{root.name}/")

    def walk(dir_path: Path, prefix=""):
        dirs = []
        files = []
        for child in sorted(dir_path.iterdir(), key=lambda p: (p.is_file(), p.name.lower())):
            if child.is_dir():
                dirs.append(child)
            else:
                files.append(child)

        for i, d in enumerate(dirs):
            is_last = (i == len(dirs) - 1) and not files
            connector = "└─ " if is_last else "├─ "
            lines.append(prefix + connector + d.name + "/")
            walk(d, prefix + ("   " if is_last else "│  "))

        for i, f in enumerate(files):
            is_last = (i == len(files) - 1)
            connector = "└─ " if is_last else "├─ "
            lines.append(prefix + connector + f.name)

    walk(root, "├─ ")
    lines.append("```")
    return "\n".join(lines)


if __name__ == "__main__":
    md_text = make_md_tree_from_fs(ROOT_DIR)

    # 出力ディレクトリが存在しない場合は作成
    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)

    # 書き込み
    OUT_PATH.write_text(md_text, encoding="utf-8")
    print(f"[OK] Markdown tree written to: {OUT_PATH}")
