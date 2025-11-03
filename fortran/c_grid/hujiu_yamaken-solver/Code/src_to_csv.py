from pathlib import Path
import csv
import os

# ======== 設定 ========
# ルートディレクトリ（走査対象）
root = Path("/root/wsl/fortran/c_grid/hujiu_yamaken-solver/Code/grid")

# 出力ディレクトリ（自由に変更可能）
out_dir = Path("/root/wsl/fortran/c_grid/hujiu_yamaken-solver/Code")  # ← 保存先を指定
# =======================

# 出力ファイル名
out_dir.mkdir(parents=True, exist_ok=True)
out_path = out_dir / "grid_inventory_with_size.csv"

rows = []
for p in root.rglob("*"):
    kind = "dir" if p.is_dir() else "file"
    level = len(p.relative_to(root).parts)
    parent = "" if p == root else p.parent.name
    ext = "" if kind == "dir" else p.suffix.lstrip(".")
    top_level = p.relative_to(root).parts[0] if level >= 1 else p.name

    # ファイルサイズとFortran行数（初期値）
    size_bytes = None
    fortran_lines = None

    if kind == "file":
        try:
            size_bytes = p.stat().st_size
            if ext in ["f90", "f", "f90inc"]:
                with open(p, "r", encoding="utf-8", errors="ignore") as f:
                    fortran_lines = sum(1 for _ in f)
        except Exception as e:
            print(f"[WARN] {p} の情報取得に失敗: {e}")

    rows.append({
        "name": p.name,
        "type": kind,
        "level": level,
        "parent": parent,
        "path": str(p.as_posix()),
        "ext": ext,
        "top_level": top_level,
        "size_bytes": size_bytes,
        "fortran_lines": fortran_lines,
    })

# ======== CSV出力 ========
with open(out_path, "w", newline="", encoding="utf-8-sig") as f:
    writer = csv.DictWriter(f, fieldnames=rows[0].keys())
    writer.writeheader()
    writer.writerows(rows)

print(f"[OK] CSVを出力しました: {out_path}")
