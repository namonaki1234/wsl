import io
import sys
from collections import defaultdict,deque,Counter

# 下記に標準入力を記載
_INPUT = """\
5 3
1 2
1 3
4 5
"""
sys.stdin = io.StringIO(_INPUT)
# ここからコードを記載
# 1-indexでやる
n,m = map(int,input().split())

gragh = [[] for _ in range(n+1)]

#隣接リスト作成
for _ in range(m):
    a,b = map(int,input().split())
    gragh[a].append(b)
    gragh[b].append(a)

# print(gragh)
visited = [False]*(n+1)

def dfs(c):
    visited[c] = True
    #cに隣接する頂点dを次に見る
    for d in gragh[c]:
        if visited[d]:
            continue
        dfs(d)

ans = 0
#すべての頂点それぞれを起点としたdfsによる全探索を行う
for i in range(1,n+1):
    if not visited[i]:
        ans += 1
        dfs(i)

print(ans)