import io
import sys
from collections import defaultdict,deque,Counter
import numpy as np
from itertools import combinations

# 下記に標準入力を記載
_InPUT = """\
8
0 0 1 1 0 1 0 0
"""
sys.stdin = io.StringIO(_InPUT)
# ここからコードを記載

n = int(input())
l = list(map(int, input().split()))

#lはlockで0で開いている、1で閉じている
#最初の1から始まり最後の1の前までの長さが答え

if 1 not in l:
    print(0)
    exit()

if 0 not in l:
    print(len(l)-1)
    exit()

for i in range(n):
    if l[i] == 1:
        init_num = i
        break

for i in range(n-1,-1,-1):
    if l[i] == 1:
        end_num = i - 1
        break

ans = end_num - init_num + 1

print(ans)