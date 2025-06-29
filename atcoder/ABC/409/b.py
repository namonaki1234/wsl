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

# xが0からnまでの値を考え、各xに対してaの要素がx以上であるものの数を数える。
for x in range(n, -1, -1):
    cnt = 0
    for a_i in a:
        if a_i >= x:
            cnt += 1
    if cnt >= x:
        print(x)
        exit()