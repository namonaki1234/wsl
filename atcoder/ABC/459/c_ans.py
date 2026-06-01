import io
import sys
from collections import Counter, deque
import itertools
import bisect
import re
# from atcoder.dsu import DSU

# 再帰の上限を増やす（Pythonでは必須のおまじない）
sys.setrecursionlimit(10**6)

# 下記に標準入力を記載
_InPUT = """\
3 7
1 1
1 3
1 3
2 1
2 2
1 2
2 1
"""
sys.stdin = io.StringIO(_InPUT)

# ここからコードを記載

n,q = map(int,input().split())
blocks = [0]*(n+1)
# height_count[j]: 削除しない世界で，j個以上あるマスの数
height_count = [0] * (q + 2)
cnt = 0
minus_cnt = 0
for i in range(q):
    cnt = 0
    query = list(map(int,input().split()))
    if query[0] == 1:
        # print(blocks)
        x = query[-1]
        blocks[x] += 1
        # blocks[x] 個以上あるマスが1つ増える
        height_count[blocks[x]] += 1
        # 全マスが blocks[x] 個以上になったなら，
        # そこまで全体削除が進んだと考えられる
        if height_count[blocks[x]] == n:
            minus_cnt = blocks[x]
    elif query[0] == 2:
        # print(blocks)
        y = query[-1]

        # 実際に y 個以上ある
        # ⇔ 削除しない世界で y + minus_cnt 個以上ある
        target = y + minus_cnt

        if target > q:
            print(0)
        else:
            print(height_count[target])

"""
なぜ q と比較するのか

クエリは全部で q 回しかない。

そのうち，ブロックを置く操作 1 x も最大で q 回である。

つまり，削除しない世界で考えても，ある1つのマスに積まれるブロック数の最大値は，

最大でも q 個

である。
"""