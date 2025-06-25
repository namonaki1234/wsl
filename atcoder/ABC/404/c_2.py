from atcoder.dsu import DSU
import io
import sys
from collections import defaultdict,deque,Counter
import math

# 下記に標準入力を記載
_INPUT = """\
4 4
2 4
3 1
4 1
2 3

"""
sys.stdin = io.StringIO(_INPUT)
# ここからコードを記載
# 入力
n, m = map(int, input().split())
deg = [0] * n
dsu = DSU(n)

for _ in range(m):
    a, b = map(int, input().split())
    a -= 1  # 0-indexed にする
    b -= 1
    deg[a] += 1
    deg[b] += 1
    dsu.merge(a, b)
    # print(f"Union: {a+1}, {b+1} -> {dsu.groups()},len{len(dsu.groups())}")  # デバッグ用出力
    # print(dsu.parent_or_size)  # デバッグ用出力
    # print(dsu.merge(a, b))  # デバッグ用出力
    # print(dsu.leader(a))  # デバッグ用出力
    # print(dsu.size(a))  # デバッグ用出力
    # print(dsu.same(a, b))  # デバッグ用出力
    print(f"Degrees: {deg}")  # デバッグ用出力

# 条件1: M == N
if m != n:
    print("No")
    exit()

# 条件2: 連結成分が1つか
if len(dsu.groups()) != 1:
    print("No")
    exit()

# 条件3: 全頂点の次数が2か
if all(d == 2 for d in deg):
    print("Yes")
else:
    print("No")
