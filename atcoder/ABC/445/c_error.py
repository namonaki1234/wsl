import io
import sys

# 下記に標準入力を記載
_InPUT = """\
7
2 4 7 5 5 6 7
"""
sys.stdin = io.StringIO(_InPUT)
# ここからコードを記載
# 最後のマス番号または最後の移動の際のx(a)を出力
# ex1_ans:5 5 7 5 5 6 7
n = int(input())
a = list(map(int,input().split()))

# s = n
# now_index = 1
# ans = []
# for i in range(1,s+1):
#     # now_index = 1
#     now_index = i
#     for j in range(1,s+1):
#         # if j == now_index:
#         #     continue
#         next_index = a[now_index-1]
#         now_index = next_index
#     ans.append(now_index)

# print(*ans)

# visited[i] : インデックス i (1..n) をすでに訪れたかどうか
# 0 番は使わないので n+1 確保する
visited = [False] * (n + 1)

# 通った順番（1-index のまま記録）
order = []

def dfs(i: int) -> None:
    """
    i : 現在地のインデックス（1-index）
    - すでに訪問済みなら終了
    - 未訪問なら記録して、a[i] に移動して再帰
    """
    # 範囲外（念のため）
    if not (1 <= i <= n):
        return

    # すでに訪問済みならここで止める（ループ防止）
    if visited[i]:
        return

    # 訪問処理
    visited[i] = True
    order.append(i)

    # 次に移動するインデックス（1-index）
    nxt = a[i - 1]      

    # 次へ（範囲内なら進む）
    if 1 <= nxt <= n:
        dfs(nxt)

# print("visit order (1-index):", order[-1])
ans = []
for i in range(1,n+1):
    # どこからスタートするか（1-index）
    start = i
    dfs(start)
    ans.append(order[-1])

print(*ans)
    