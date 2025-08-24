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
s = [input().strip() for _ in range(n)]

#pt=pointのこと
pt = [0] * n

for j in range(m):
    x = sum(row[j] == '0' for row in s)  # 0 の数
    y = n - x                            # 1 の数
    for i in range(n):
        #N が奇数であることからx=y となることはない、かつx=0 またはy=0 である場合、全員に1 点が与えられるので、 x <= y:と x >= y:のように両方に=をつけて問題ない
        if s[i][j] == '0':
            if x <= y:
                pt[i] += 1
        else:
            if x >= y:
                pt[i] += 1

high = max(pt)
ans = [str(i + 1) for i in range(n) if pt[i] == high]
print(" ".join(ans))