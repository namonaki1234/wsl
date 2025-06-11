import io
import sys
from collections import defaultdict,deque,Counter
import numpy as np

# 下記に標準入力を記載
_InPUT = """\
12
11 3 8 9 5 2 10 4 1 6 12 7

"""
sys.stdin = io.StringIO(_InPUT)
# ここからコードを記載

n = int(input())
p = list(map(int, input().split()))
copy_p = p.copy()

for i in range(n-1):
    if p[i] < p[i+1] and i < n-1:
        copy_p[i]="<"
    elif p[i] > p[i+1]:
        copy_p[i]=">"


print("".join(map(str, copy_p[:-1])))


