import io
import sys
from collections import defaultdict,deque,Counter
import math

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
N, M = map(int, input().split())
edges = [set(map(int, input().split())) for _ in range(N)]

def is_single_cycle_graph(N, edges):
    if len(edges) != N:
        return False  # 条件2: 辺数が N でないと閉路にはならない

    graph = defaultdict(set)
    for u, v in edges:
        graph[u].add(v)
        graph[v].add(u)

    # 条件3: 各頂点の次数が2かどうか
    for v in graph:
        if len(graph[v]) != 2:
            return False

    # 条件1: 全頂点が連結かどうか → BFS or DFS
    visited = set()
    def bfs(start):
        queue = deque([start])
        while queue:
            v = queue.popleft()
            if v in visited:
                continue
            visited.add(v)
            for neighbor in graph[v]:
                if neighbor not in visited:
                    queue.append(neighbor)

    bfs(start=next(iter(graph)))  # 任意の1ノードからスタート
    return len(visited) == N

if (is_single_cycle_graph(N, edges)):
    print("Yes")
else:
    print("No")
