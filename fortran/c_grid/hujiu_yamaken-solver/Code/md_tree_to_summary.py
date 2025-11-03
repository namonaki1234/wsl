"""
md_tree_to_summary.py
Markdown形式のディレクトリツリー（例: treeコマンド出力や .mdの```textブロック）を読み込み、
拡張子ごとの件数とファイルパス一覧をCSVとDataFrameで出力する。
"""

import re
import pandas as pd
from pathlib import PurePosixPath, Path
from collections import defaultdict

# ======= 設定セクション =======
MD_PATH = Path("/root/wsl/fortran/c_grid/hujiu_yamaken-solver/Code/solver_tree.md")       # ← 解析対象の .md ファイル
OUT_CSV = Path("/root/wsl/fortran/c_grid/hujiu_yamaken-solver/Code/solver_extension_summary.csv")  # ← 出力先
# =================================

def parse_tree_markdown(md_text: str):
    """Markdownのツリー表記 (├─, └─, │) をパースしてファイルパスのリストを返す"""
    content = re.sub(r"^```.*\n|\n```$", "", md_text.strip(), flags=re.MULTILINE)
    lines = [ln.rstrip() for ln in content.splitlines() if ln.strip()]

    stack = []
    file_paths = []

    # 1行目（root）を処理
    if lines and lines[0].endswith("/"):
        root = lines[0].rstrip("/")
        stack = [root]
        start_idx = 1
    else:
        root = "root"
        stack = [root]
        start_idx = 0

    pat = re.compile(r"^(?P<prefix>[\s│]*)(?:├─ |└─ )(?P<name>.+)$")

    for raw in lines[start_idx:]:
        m = pat.match(raw)
        if not m:
            fallback = re.sub(r"^[\s│├─└]+", "", raw).strip()
            if not fallback:
                continue
            name = fallback
        else:
            name = m.group("name").strip()

        if name.endswith("/"):
            stack.append(name.rstrip("/"))
        else:
            # フルパス生成
            path = PurePosixPath("/".join(stack)) / name
            file_paths.append(path)

    return file_paths


def get_ext(p: PurePosixPath):
    sfx = p.suffix.lower()
    return sfx if sfx else "[noext]"


def make_summary(md_path: Path):
    text = md_path.read_text(encoding="utf-8", errors="ignore")
    paths = parse_tree_markdown(text)

    ext_to_paths = defaultdict(list)
    for p in paths:
        ext_to_paths[get_ext(p)].append(str(p))

    rows = []
    for ext, plist in ext_to_paths.items():
        rows.append({
            "extension": ext,
            "count": len(plist),
            "paths": "\n".join(sorted(plist, key=lambda x: x.lower()))
        })

    df = pd.DataFrame(rows).sort_values(["count", "extension"], ascending=[False, True])
    return df.reset_index(drop=True)


if __name__ == "__main__":
    df = make_summary(MD_PATH)
    OUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUT_CSV, index=False, encoding="utf-8")
    print(f"[OK] 拡張子ごとの集計結果を出力しました: {OUT_CSV}")
    print(df)
