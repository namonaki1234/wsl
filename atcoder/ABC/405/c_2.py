import io
import sys
from collections import defaultdict,deque,Counter
import numpy as np
import itertools

# 下記に標準入力を記載
_InPUT = """\
10
7781 8803 8630 9065 8831 9182 8593 7660 7548 8617

"""
sys.stdin = io.StringIO(_InPUT)
# ここからコードを記載

n = int(input())
a = list(map(int, input().split()))

sum_to_j_1 = 0
ans = 0
for i in range(n):
    if n-1 == i:
        break
    sum_to_j_1 += a[i]
    ans += a[i+1] * sum_to_j_1


# print(a_sumproducts)
# print(sum(a_sumproducts))
print(ans)