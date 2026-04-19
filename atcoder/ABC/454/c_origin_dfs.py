import io
import sys
import itertools
from atcoder.dsu import DSU
# 再帰の上限を増やす（Pythonでは必須のおまじない）
sys.setrecursionlimit(10**6)

# 下記に標準入力を記載
_InPUT = """\
5 5
1 2
2 3
3 4
2 4
5 2
"""
sys.stdin = io.StringIO(_InPUT)
# ここからコードを記載
#n種類のアイテム、はじめ高橋君はアイテム1を持っている
#高橋にはM人と友達がいる。i人目の友達にアイテムA_iを渡すと、アイテムB_iをもらうことができる
#高橋が手に入れることのできるアイテムはアイテム1を含めて何種類

#今回は1-indexだから入力変数もn+1で、初期値も1
#もし0-indexなら入力変数をn、初期値を0としなくてはいけない
n,m = map(int, input().split())

ab = []
for _ in range(m):
    ab.append([int(x) for x in input().split()])

# print(ab)

# takahashi_item = [1]

def dfs(current_item):
    for next_item in gragh[current_item]:
        if not visited[next_item]:
            visited[next_item] = True
            dfs(next_item)

gragh = [[] for _ in range(n+1)]
edges = ab
for a,b in edges:
    gragh[a].append(b)

visited = [False] * (n+1)
visited[1] = True

dfs(1)

print(sum(visited))