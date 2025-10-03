import io
import sys

# 下記に標準入力を記載 今回のエラーの原因がわかるような例を考えた
_InPUT = """\
6
6 1
2 3
3 1
0 0
4 6
0 0
"""
sys.stdin = io.StringIO(_InPUT)
# ここからコードを記載

from collections import deque

N = int(input())
G = [[] for _ in range(N+1)]

# グラフを構築
for i in range(1, N+1):
    a, b = map(int, input().split())
    G[a].append(i)
    G[b].append(i)

ok = [0] * (N+1)
ok[0] = 1  # スキル0は最初から習得済み

# BFS開始
q = deque([0])
while q:
    v = q.popleft()
    for vv in G[v]:
        if not ok[vv]:
            ok[vv] = 1
            q.append(vv)

print(sum(ok) - 1)  # スキル0を除いてカウント

