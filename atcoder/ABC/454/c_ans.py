import io
import sys
from collections import deque
import itertools
from atcoder.dsu import DSU
# 再帰の上限を増やす（Pythonでは必須のおまじない）
sys.setrecursionlimit(10**6)

# 下記に標準入力を記載
_InPUT = """\
3 2
2 1
3 2
"""
sys.stdin = io.StringIO(_InPUT)
# ここからコードを記載
#n種類のアイテム、はじめ高橋君はアイテム1を持っている
#高橋にはM人と友達がいる。i人目の友達にアイテムA_iを渡すと、アイテムB_iをもらうことができる
#高橋が手に入れることのできるアイテムはアイテム1を含めて何種類

n,m = map(int, input().split())
g = [[] for _ in range(n)]
for _ in range(m):
    x, y = map(int, input().split())
    g[x - 1].append(y - 1)
q = deque()      # 予約リスト（次にどのアイテムを起点に交換を試みるか）
d = [False] * n  # 記憶（どのアイテムをすでに手に入れたか）

q.append(0)      # 最初はアイテム1（インデックス0）を持っている
d[0] = True      # アイテム1は手に入れたと記憶する
while q:
    x = q.popleft()    # 待合室からアイテムを1つ取り出す
    for i in g[x]:     # そのアイテム(x)から交換できるアイテム(i)をすべて確認
        if d[i]:       # すでに持っているアイテムならスルー
            continue
        d[i] = True    # 新しく手に入れた！とノートに記録
        q.append(i)    # この新しいアイテム(i)を使って、さらに別の交換ができるか後で調べる
print(d.count(True))
