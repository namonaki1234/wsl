import io
import sys
from collections import defaultdict,deque,Counter
import numpy as np
from itertools import combinations

# 下記に標準入力を記載
_InPUT = """\
5 5
2 3
1 2 1000000
3 4
2 2
2 3

"""
sys.stdin = io.StringIO(_InPUT)
# ここからコードを記載

n, q = map(int, input().split())
query = [list(map(int, input().split())) for _ in range(q)]

a = [i for i in range(1, n + 1)]

offset = 0
for i in range(q):
    if len(query[i]) != 2:
        t,p,x = query[i]
        a[(p-1+offset)%n] = x
    else:
        t, p_or_k = query[i]
        if t == 2:
            print(a[(p_or_k-1+offset)%n])
        elif t == 3:
            offset += p_or_k%n
