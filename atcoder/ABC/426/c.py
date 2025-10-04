import io
import sys
from collections import defaultdict,deque,Counter

# 下記に標準入力を記載
_InPUT = """\
8 5
2 6
3 5
1 7
5 7
7 8

"""
sys.stdin = io.StringIO(_InPUT)
# ここからコードを記載
sys.setrecursionlimit(10**9)

n, q = map(int,input().split())

# os = [i for i in range(1,n+1)]
# print(os)

G=[i for i in range(1,n+1)]
# print(G)
cnt = 0
for i in range(1,q+1):
  cnt = 0
  x,y=map(int,input().split())
  for j in range(1,x+1):
    if G[j-1] <= x:
      G[j-1] = y
      cnt += 1
#   G[x-1].append(i)
#   G[y-1].append(i)
  print(cnt)
# print(G)
# ok=[0]*(q+1)
# ok[0]=1

# def dfs(v):
#   ok[v]=1
#   for vv in G[v]:
#     if not ok[vv]:
#       dfs(vv)

# dfs(0)

# print(sum(ok)-1)





# cnt = 0
# for _ in range(q):
#     cnt = 0
#     x,y = map(int,input().split())

#     os.sort()
#     for i in os:
#         if i > x:
#             break
#         if i <= x:
#             os[os.index(i)] = y
#             cnt += 1
#     print(os)
#     print(cnt)

# cnt = 0
# for _ in range(q):
#     cnt = 0
#     x,y = map(int,input().split())
#     for i in range(len(os)):
#         if os[i] <= x:
#             os[i] = y
#             cnt += 1
#     print(cnt)
