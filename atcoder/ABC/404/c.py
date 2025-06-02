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
graph = defaultdict(set)
edges = [set(map(int, input().split())) for _ in range(M)]

for u, v in edges:
    graph[u].add(v)
    graph[v].add(u)

# print(graph)

Counter_graph = Counter(graph)
# print(Counter_graph)

if N != M:
    print("No")
    exit()

for k, v in Counter_graph.items():
    # print(k, v)
    if len(v) != 2:
        print("No")
        exit()

print("Yes")
