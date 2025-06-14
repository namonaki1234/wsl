import io
import sys
from collections import defaultdict,deque,Counter
import numpy as np
from itertools import combinations

# 下記に標準入力を記載
_InPUT = """\
10 12
4 4 5 7 1 7 0 8 5
"""
sys.stdin = io.StringIO(_InPUT)
# ここからコードを記載

n, l = map(int, input().split())
d = list(map(int, input().split()))

p = [0]
for i in range(n-1):
    p.append((p[i] + d[i])%6)

com = list(combinations(p, 3))
com_sorted = [tuple(sorted(t)) for t in com]
print(com_sorted)
print(p)

count = 0
if l% 3 != 0:
    print(0)
    exit()
for i in range(len(com)):
    # a, b, c = com[i]
    # if d[b-1] - d[a-1] == d[c-1]- d[b-1]:
    #     count += 1
    #     print(count)
    # if (com_sorted[i][1] - com_sorted[i][0]) == l/3 and (com_sorted[i][2] - com_sorted[i][1]) == l/3 and sum(com_sorted[i]) == l:
    if (com_sorted[i][1] - com_sorted[i][0]) == l/3 and (com_sorted[i][2] - com_sorted[i][1]) == l/3 and (com_sorted[i][2] - com_sorted[i][0]) == 2*l/3:
        count += 1

print(count)