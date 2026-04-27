import io
import sys
from collections import Counter, deque
import itertools
# from atcoder.dsu import DSU

# 再帰の上限を増やす（Pythonでは必須のおまじない）
sys.setrecursionlimit(10**6)

# 下記に標準入力を記載
_InPUT = """\
10 2
3 3 4 1 1 3 3 1 5 1
"""
sys.stdin = io.StringIO(_InPUT)

# ここからコードを記載

n, k = map(int, input().split())
a = [int(x) for x in input().split()]

# 各値が何回出てくるか数える
counter = Counter(a)

# 最初の合計値
total = sum(a)

# その値を選んだときに消える合計値を入れる
vanish_values = []

for x, count in counter.items():
    vanish_values.append(x * count)

# 大きい順に並べる
vanish_values.sort(reverse=True)

# 大きいものから最大 k 個消す
for v in vanish_values[:k]:
    total -= v

print(total)