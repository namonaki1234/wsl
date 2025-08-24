import io
import sys
from collections import defaultdict,deque,Counter
import numpy as np
from itertools import combinations

# 下記に標準入力を記載
_InPUT = """\
3 5
11100
10101
01110

"""
sys.stdin = io.StringIO(_InPUT)
# ここからコードを記載

n,m = map(int, input().split())

for r in range(n):
    s = input().strip()
    s_bit = int(s, 2)
    print(s_bit)
    if s_bit == 0:
        print(-1)
    else:
        index = m - s_bit.bit_length()
        print(r, index)