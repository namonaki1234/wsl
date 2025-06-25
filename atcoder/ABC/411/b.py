import io
import sys
from collections import defaultdict,deque,Counter
import numpy as np

# 下記に標準入力を記載
_InPUT = """\
2
100

"""
sys.stdin = io.StringIO(_InPUT)
# ここからコードを記載

n = int(input())
d = list(map(int, input().split()))

sum =[0] * (n - 1)
sum[0] = d[0]
for i in range(1,n-1):
    sum[i] += sum[i-1] + d[i]
# print(*d)

print(*sum)
temp =[]
for i in range(0,n-1):
    for j in range(i+1, n-1):
        temp.append(sum[j] - sum[i])
    print(*temp)
    temp = []
