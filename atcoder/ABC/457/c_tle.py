import io
import sys
from collections import Counter, deque
import itertools
import bisect
# from atcoder.dsu import DSU

# 再帰の上限を増やす（Pythonでは必須のおまじない）
sys.setrecursionlimit(10**6)

# 下記に標準入力を記載
_InPUT = """\
3 9
3 1 3 2
1 3
2 4 3
1 3 2
"""
sys.stdin = io.StringIO(_InPUT)

# ここからコードを記載

n,k = map(int,input().split())

a = []
for i in range(n):
    row = [z for z in input().split()][1:]
    a.append(row)

c = [y for y in input().split()]

b = []
for i in range(len(c)):
    for j in range(int(c[i])):
        b.append(a[i])

# print(b)

# 各b配列つまりrowの開始位置をmemo (累積和)
prefix_sum = [0]
for row in b:
    prefix_sum.append(prefix_sum[-1] + len(row))

target_idx = k - 1

# prefix_sum内でtargetが含まれるrowが何番目(idx+1)かを示す
row_idx = bisect.bisect_right(prefix_sum, target_idx) - 1

# prefix_sum[row_idx]はtargetが含まれるrowの開始idx
col_idx = target_idx - prefix_sum[row_idx]

# print(row_idx,col_idx,b[row_idx][col_idx])
print(b[row_idx][col_idx])