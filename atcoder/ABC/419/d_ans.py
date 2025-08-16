import io
import sys
from collections import defaultdict,deque,Counter
import numpy as np
from itertools import combinations,permutations

# 下記に標準入力を記載
_InPUT = """\
5 3
apple
lemon
2 4
1 5
5 5


"""
sys.stdin = io.StringIO(_InPUT)
# ここからコードを記載
# sys.setrecursionlimit(1000000)

n,m = map(int, input().split())
s = list(input())
t = list(input())

lr_imos_cnt = [0] * (n + 1)
for _ in range(m):
    l, r = map(int, input().split())
    lr_imos_cnt[l-1] += 1
    lr_imos_cnt[r] -= 1

# # 累積和の準備
# for i in range(n):
#     lr_imos_cnt[i+1] += lr_imos_cnt[i]
#     print(lr_imos_cnt[i])

# print(lr_imos_cnt)

# 累積和
for i in range(n):
    lr_imos_cnt[i+1] += lr_imos_cnt[i]
    # print(lr_imos_cnt[i])

# print(lr_imos_cnt)

for i in range(n):
    if lr_imos_cnt[i] % 2 == 0:
        # print(i, lr_imos_cnt[i])
        continue
    s[i], t[i] = t[i], s[i]
print("".join(s))