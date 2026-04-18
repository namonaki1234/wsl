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


#もしゼロインデックスならnにして、初期値を0にすればいい
n,m = map(int, input().split())

gragh = [[] for _ in range(n+1)]

ab = []
for _ in range(m):
    ab.append([int(x) for x in input().split()])

edges = ab
for a,b in edges:
    gragh[a].append(b)

q = deque()
q.append(1)
d = [False] * (n+1)
d[1] = True

while q:
    x = q.popleft()
    for i in gragh[x]:
        if d[i]:
            continue
        d[i] = True
        q.append(i)

print(d.count(True))