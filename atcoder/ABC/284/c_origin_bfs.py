import io
import sys
from collections import defaultdict,deque,Counter

# 下記に標準入力を記載
_INPUT = """\
4 6
1 2
1 3
1 4
2 3
2 4
3 4
"""
sys.stdin = io.StringIO(_INPUT)
# ここからコードを記載
# 1-indexでやる
n,m = map(int,input().split())

gragh = [[] for _ in range(n+1)]

for _ in range(m):
    a,b = map(int,input().split())
    gragh[a].append(b)
    gragh[b].append(a)

visited = [False]*(n+1)
ans = 0
q = deque()

#すべての頂点それぞれについてbfsによる全探索
for i in range(1,n+1):
    if visited[i]:
        continue

    ans += 1
    visited[i] = True
    #初期値設定
    q.append(i)

    while q:
        c = q.popleft()
        for d in gragh[c]:
            if visited[d]:
                continue
            visited[d] = True
            q.append(d)

print(ans)