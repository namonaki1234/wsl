import io
import sys
from collections import Counter, deque
import itertools
# from atcoder.dsu import DSU

# 再帰の上限を増やす（Pythonでは必須のおまじない）
sys.setrecursionlimit(10**6)

# 下記に標準入力を記載
_InPUT = """\
8 6
1 2 3 4 1 2 3 4
"""
sys.stdin = io.StringIO(_InPUT)

# ここからコードを記載

n, k = map(int, input().split())
a = [int(x) for x in input().split()]

# 各値が何回出てくるか数える
counter_a = Counter(a)

delete_sum = []
for key,value in counter_a.items():
    delete_sum.append(key * value)

delete_sum.sort(reverse=True)
# print(delete_sum)
total = sum(a)

if k > len(delete_sum):
    k = len(delete_sum)

for i in range(k):
    total -= delete_sum[i]
print(total)
