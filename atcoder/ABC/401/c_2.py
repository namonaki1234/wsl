import io
import sys
from collections import defaultdict,deque,Counter
import math
# from more_itertools import chunked

# 下記に標準入力を記載
_INPUT = """\
4 2

"""
sys.stdin = io.StringIO(_INPUT)
# ここからコードを記載
n, k = map(int, input().split())

a =[1 for i in range(n + 1)]  # a[i]はi番目の数を表す
# 0からk-1までの数の和がkとなる。なぜなら、a[0]からa[k-1]までの数はすべて1である。
s = k

for i in range(k, n + 1):
    a[i] = s
    s -= a[i - k]
    s += a[i]
    s %= 1000000000  # 1000000000で割った余りを求める

print(a[n])