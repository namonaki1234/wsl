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
blocks = [0]*n
cnt = 0
minus_cnt = 0
for i in range(q):
    cnt = 0
    query = list(map(int,input().split()))
    if query[0] == 1:
        # print(blocks)
        x = query[-1]
        blocks[x-1] += 1
        if min(blocks) >=1:
            for j in range(n):
                # blocks[j] -=1
                minus_cnt += 1
    elif query[0] == 2:
        # print(blocks)
        y = query[-1]
        # print("test",min(blocks),y)
        for k in range(n):
            if (blocks[k]-minus_cnt) >= y:
                cnt += 1
                # print(blocks[k])
                # print(blocks)
        print(cnt)