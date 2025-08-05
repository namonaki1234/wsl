import io
import sys
from collections import defaultdict,deque,Counter
import numpy as np
from itertools import combinations,permutations
import math

# 下記に標準入力を記載
_InPUT = """\
9
3 1 4 1 5 9 2 6 5

"""
sys.stdin = io.StringIO(_InPUT)
# ここからコードを記載
sys.setrecursionlimit(1000000)

n = int(input())
a = list(map(int, input().split()))

ans = 0
# ans = set()
# visited = [[False] * (n+1) for _ in range(n+1)]
visited = [False] * (n + 1)
def dfs(i,j):
    global ans
    # visited[i][j] = True
    visited[i] = True
    if i >= j or j > n:
        return
    # if j - i == a[i-1]+a[j-1] and not visited[j][i]:
    if j - i == a[i-1]+a[j-1] and not visited[j]:
        ans += 1
        # ans.add((i,j))
        # print(i,j,a[i-1]+a[j-1])
    print(i,j,a[i-1]+a[j-1])
    # for _ in range(1,n+1):
    #     dfs(i,j-1)
    # dfs(i,j+1)
    dfs(i,j-1)
    dfs(i+1,j)
dfs(1,n)

# print(len(ans))
print(ans)