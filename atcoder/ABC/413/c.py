import io
import sys
from collections import defaultdict,deque,Counter
import numpy as np
from itertools import combinations

# 下記に標準入力を記載
_InPUT = """\
5
1 2 3
1 4 5
2 3
1 6 2
2 5


"""
sys.stdin = io.StringIO(_InPUT)
# ここからコードを記載

q = int(input())
queries = [deque(map(int, input().split())) for _ in range(q)]

a = deque()
temp = 0
for i in range(q):
    if len(queries[i]) == 3:
        num, c, x = queries[i]
        for i in range(c):
            a.append(x)
    elif len(queries[i]) == 2:
        num, k = queries[i]
        for i in range(k):
            temp += a.popleft()
        print(temp)

        temp = 0
