import io
import sys
from collections import defaultdict,deque,Counter
import numpy as np
import itertools as itbbbE

# 下記に標準入力を記載
_InPUT = """\
1 2
1
1 1


"""
sys.stdin = io.StringIO(_InPUT)
# ここからコードを記載

n, m = map(int, input().split())
a = list(map(int, input().split()))
b = list(map(int, input().split()))

ans = a.copy()
for item_a in a:
    for item_b in b:
        if item_a == item_b:
            ans.remove(item_b)
            b.remove(item_b)
            break

# def dfs(n, use3, use5, use7):
#     global ans
#     if n > N:
#         return
#     if use3 and use5 and use7:
#         ans += 1
#     dfs(10*n+3, True, use5, use7)
#     dfs(10*n+5, use3, True, use7)
#     dfs(10*n+7, use3, use5, True)

# dfs(0, False, False, False)

print(*ans)