import io
import sys
from collections import defaultdict,deque,Counter

# 下記に標準入力を記載
_INPUT = """\
5
1 futuremap
1 howtospeak
2
3
2
"""
sys.stdin = io.StringIO(_INPUT)


Q = int(input())
queries = [input().split() for _ in range(Q)]
# print(queries)

s = deque()
for q in queries:
    if q[0] == "1":
        s.append(q[-1])
    elif q[0] == "2":
        print(s[-1])
    elif q[0] == "3":
        s.pop()