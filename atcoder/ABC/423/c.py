import io
import sys
from collections import defaultdict,deque,Counter
import numpy as np
from itertools import combinations

# 下記に標準入力を記載
_InPUT = """\
6 3
1 0 0 1 0 0
"""
sys.stdin = io.StringIO(_InPUT)
# ここからコードを記載

n,r = map(int, input().split())
l = list(map(int, input().split()))

#lはlockで0で開いている、1で閉じている
#rの場所から始まって、すべて1にするまでの最小回数を求める

if 1 not in l:
    print(len(l))
    exit()

if 0 not in l:
    print(0)
    exit()

for i in range(n):
    if l[i] == 0:
        init_num = i
        break

for i in range(n-1,-1,-1):
    if l[i] == 0:
        end_num = i
        break

ans = end_num - init_num + min(abs(r-end_num),abs(r-init_num))


print(ans)