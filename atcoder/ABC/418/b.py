import io
import sys
from collections import defaultdict,deque,Counter
import numpy as np
from itertools import combinations

# 下記に標準入力を記載
_InPUT = """\
ottottott


"""
sys.stdin = io.StringIO(_InPUT)
# ここからコードを記載

s = list(input().strip())
t_indices = [i for i,j in enumerate(s) if j == 't']

# print(t_indices)
cnt = 0
ans = 0

t_comb = list(combinations(t_indices, 2))

# print(t_comb)

for i, j in t_comb:
    if j - i + 1 >= 3:
        s_slice = s[i:j+1]
        x = s_slice.count('t')
        ans = max(ans, (x-2)/(j-i+1-2))
        # print(i, j, x, ans)

print(ans)
