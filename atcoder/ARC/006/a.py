import io
import sys
from collections import defaultdict,deque,Counter
import numpy as np

# 下記に標準入力を記載
_InPUT = """\
1 2 3 4 5 6
7
1 2 3 4 5 6
"""
sys.stdin = io.StringIO(_InPUT)
# ここからコードを記載

e = list(map(int, input().split()))
b = int(input())
l = list(map(int, input().split()))

cnt = 0
for e_i in e:
    for l_i in l:
        if e_i == l_i:
            cnt += 1

if cnt == 6:
    print("1")
elif cnt == 5:
    print("2")
elif cnt == 4:
    print("3")
elif cnt == 3:
    print("4")
