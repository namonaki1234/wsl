import io
import sys
from collections import defaultdict,deque,Counter
import numpy as np

# 下記に標準入力を記載
_InPUT = """\
7
1 6 2 10 2 3 2
"""
sys.stdin = io.StringIO(_InPUT)
# ここからコードを記載

n = int(input())
a = list(map(int, input().split()))

a_sorted = sorted(a, reverse=True)
# print(a_sorted)

ans = []
for i in range(n):
    # if a_sorted[i] <= i+1:
    #     ans = a_sorted[i]
    #     print(ans)
    #     exit()
# if len(ans) == 0:
#     ans = 0
#     print(ans)
#     exit()
    if (i+1)-a_sorted[i] == 0:
        ans = a_sorted[i]
        print(ans)
        exit()
    if (i+1)-a_sorted[i] == 1:
        ans = a_sorted[i]
        print(ans)
        exit()
    if (i+1)-a_sorted[i] >= 2:
        ans = i
        print(ans)
        exit()
