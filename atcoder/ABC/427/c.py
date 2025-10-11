import io
import sys
from collections import defaultdict,deque,Counter

# 下記に標準入力を記載
_InPUT = """\
5 8
1 2
1 3
1 4
2 3
2 5
3 4
3 5
4 5


"""
sys.stdin = io.StringIO(_InPUT)
# ここからコードを記載
sys.setrecursionlimit(10**9)

n, m = map(int, input().split())
edges = [tuple(map(int, input().split())) for _ in range(m)]

graph = defaultdict(list)
for a, b in edges:
    graph[a].append(b)
    graph[b].append(a)

print(graph[1])
visited = [False] * (n + 1)
ans = 0

def dfs(v):
    visited[v] = True
    for nv in graph[v]:
        if not visited[nv]:
            dfs(nv)

for i in range(1, n + 1):
    if not visited[i]:
        dfs(i)
        ans += 1

print(ans)