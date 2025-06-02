import io
import sys
from collections import deque

# 下記に標準入力を記載
_INPUT = """\
3 5
.#?#.
.?#?.
?...?


"""
sys.stdin = io.StringIO(_INPUT)
# ここからコードを記載

H, W = map(int, input().split())
grid = [list(input().strip()) for _ in range(H) ]

rows, cols = H, W

for r in range(rows):
    for c in range(cols):
        if grid[r][c] == '?':
            grid[r][c] = grid[r][c].replace("?", "#")             
        else:
            continue
 



visited = set()  # 訪問済みのセルを記録するセット

# グリッド内の最初の`#`を見つける
for r in range(rows):
    for c in range(cols):
        if grid[r][c] == '#':
            start_row, start_col = r, c
            break
    else:
        continue
    break

# BFSの初期設定
queue = deque([(start_row, start_col)])
visited.add((start_row, start_col))

min_row, max_row = start_row, start_row
min_col, max_col = start_col, start_col

# BFSを実行
while queue:
    row, col = queue.popleft()

    # 長方形の範囲を更新
    min_row = min(min_row, row)
    max_row = max(max_row, row)
    min_col = min(min_col, col)
    max_col = max(max_col, col)

    # 隣接するセルを探索
    for dr, dc in [(-1, 0), (1, 0), (0, -1), (0, 1)]:  # 上下左右
        new_row, new_col = row + dr, col + dc
        if 0 <= new_row < rows and 0 <= new_col < cols:
            if (new_row, new_col) not in visited and grid[new_row][new_col] == '#':
                visited.add((new_row, new_col))
                queue.append((new_row, new_col))

# 結果を出力
#print(f"長方形の範囲: 上 {min_row}, 下 {max_row}, 左 {min_col}, 右 {max_col}")

#square = (max_row - min_row + 1) * (max_col - min_col + 1)

sorted_visited = sorted(visited)
 
for i in range(min_row, max_row+1):
    for j in range(min_col, max_col+1):
        if grid[i][j] != "#" :
            print("No")
            exit()
            
print("Yes")




