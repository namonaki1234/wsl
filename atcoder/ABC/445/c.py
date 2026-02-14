import io
import sys

# 下記に標準入力を記載
_InPUT = """\
15
11 3 10 7 15 10 10 11 11 13 11 12 14 14 15

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

# いったん作っておく（中身はループごとに差し替える）
visited = [False] * (n + 1)
order = []

def dfs(i: int) -> None:
    if not (1 <= i <= n):
        return
    if visited[i]:
        return

    visited[i] = True
    order.append(i)

    nxt = a[i - 1]
    if 1 <= nxt <= n:
        dfs(nxt)

ans = []
for start in range(1, n + 1):
    # ★ここだけが本質的な修正：スタートごとにリセット
    visited = [False] * (n + 1)
    order = []

    dfs(start)
    ans.append(order[-1])  # この start に対するゴール

print(*ans)
    