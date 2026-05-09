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


# block_ends[i] = i番目のブロックまでの長さ
block_ends = [0]
total_sum = 0

for i in range(n):
    # 「行の長さ * 繰り返す回数」を一気に足す（ループを回さない！）
    block_total_elements = len(a[i]) * int(c[i])
    total_sum += block_total_elements
    block_ends.append(total_sum)

# print(a)
# print(block_ends)

target_idx = k - 1 # idxにする

# ターゲットが「何番目のブロック(a[i])」にあるか探す
# bisect_right でターゲットのブロックを特定
block_idx = bisect.bisect_right(block_ends, target_idx) - 1

# print(block_idx)

# そのブロックの開始位置を引いて、ブロック内での相対位置を算出
relative_idx = target_idx - block_ends[block_idx]

# print(relative_idx)

# そのブロック（a[block_idx]）の長さで割った余りが、a[block_idx] 内のidx→c[i]回繰り返しているため
ans_idx = relative_idx % len(a[block_idx])

print(a[block_idx][ans_idx])