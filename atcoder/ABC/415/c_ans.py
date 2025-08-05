import io
import sys
from collections import defaultdict,deque,Counter
import numpy as np
from itertools import combinations,permutations
import math

# 下記に標準入力を記載
_InPUT = """\
5
3
0010000
3
0010110
1
1
2
100
4
001110010101110

"""
sys.stdin = io.StringIO(_InPUT)
# ここからコードを記載

t = int(input())
for _ in range(t):
    n = int(input())
    s = input()
    s = '0' + s  # 先頭に0を追加（1-based indexにするため）

    ok = [0] * (1 << n)
    ok[0] = 1  # 初期状態

    for i in range(1 << n):
        if ok[i] == 0:
            continue
        for j in range(n):
            if i & (1 << j):
                continue
            next_state = i | (1 << j)
            if s[next_state] == '0':
                ok[next_state] = 1

#0と1はtrueとfalseになるので、それを利用して出力
    print("Yes" if ok[(1 << n) - 1] else "No")
