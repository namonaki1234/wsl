import io
import sys
from collections import defaultdict,deque,Counter
import numpy as np
from itertools import combinations

# 下記に標準入力を記載
_InPUT = """\
3 2 5
1 1
3 2
2 1
3 1
1 2

"""
sys.stdin = io.StringIO(_InPUT)
# ここからコードを記載

n,m,k = map(int,input().split())

# points = [0]*n
# ans = []
# for _ in range(k):
#     a,b = map(int,input().split())
#     points[a-1] += 1
#     if max(points) == m:
#         # ans.append(points.index(max(points))+1)
#         ans.append(points.index(points[a-1])+1)

# if max(points) != m:
#     exit()
# # print(points)
# print(*ans)
ans = []
questions = [[0]*m for _ in range(n)]
for _ in range(k):
    a,b = map(int,input().split())
    questions[a-1][b-1] += 1
    if sum(questions[a-1]) == m:
        ans.append(a)
# print(questions)
print(*ans)