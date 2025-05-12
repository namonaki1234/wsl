import io
import sys
from collections import defaultdict,deque,Counter
import numpy as np
import itertools

# 下記に標準入力を記載
_InPUT = """\
10
7781 8803 8630 8831 9182 8593 7660 7548 8617

"""
sys.stdin = io.StringIO(_InPUT)
# ここからコードを記載

n = int(input())
a = list(map(int, input().split())) 

combinations = list(itertools.combinations(range(0,len(a)), 2))

times = 0
for i in range(len(combinations)):
    a1 = combinations[i][0]
    a2 = combinations[i][1]
    times += a[a1] * a[a2]

print(times)

    


    











