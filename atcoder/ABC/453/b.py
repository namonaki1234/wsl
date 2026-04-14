import io
import sys
from collections import Counter

# 下記に標準入力を記載
_InPUT = """\
6 10
30 35 40 21 30 12 31
"""
sys.stdin = io.StringIO(_InPUT)
# ここからコードを記載

t,x = map(int, input().split())
a = [int(x) for x in input().split(" ")]

before = 0

for i in range(t+1):
    if i == 0:
        before = a[i]
        print(i,a[i])
    else:
        if (abs(a[i] - before)) >=x:
            before = a[i]
            print(i,a[i])
        # else:
        #     print(i,a[i])

