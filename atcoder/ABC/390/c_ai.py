import io
import sys
import networkx as nx

# 下記に標準入力を記載
_INPUT = """\
3 5
.#?#.
.?#?.
?...?
"""
sys.stdin = io.StringIO(_INPUT)

H, W = map(int, input().split())
grid = [list(input().strip()) for _ in range(H)]

# グラフの作成
G = nx.Graph()

# グリッドの各セルをノードとして追加
for r in range(H):
    for c in range(W):
        if grid[r][c] == '#':
            G.add_node((r, c))
            # 隣接するセルをエッジとして追加
            for dr, dc in [(-1, 0), (1, 0), (0, -1), (0, 1)]:  # 上下左右
                nr, nc = r + dr, c + dc
                if 0 <= nr < H and 0 <= nc < W and grid[nr][nc] == '#':
                    G.add_edge((r, c), (nr, nc))

# 最初の`#`を見つける
start_node = next(((r, c) for r in range(H) for c in range(W) if grid[r][c] == '#'), None)

if start_node is None:
    print("No")
    sys.exit()

# BFSを実行
visited = set(nx.bfs_tree(G, start_node))

# 長方形の範囲を計算
min_row = min(r for r, c in visited)
max_row = max(r for r, c in visited)
min_col = min(c for r, c in visited)
max_col = max(c for r, c in visited)

# 結果を出力
for i in range(min_row, max_row + 1):
    for j in range(min_col, max_col + 1):
        if grid[i][j] != "#":
            print("No")
            sys.exit()

print("Yes")
