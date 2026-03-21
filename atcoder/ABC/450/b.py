import io
import sys
from collections import Counter

# 下記に標準入力を記載
_InPUT = """\
4
25 40 65
30 55
25
"""
# 5
# 25 40 65 80
# 30 55 60
# 25 30
# 20
# """
sys.stdin = io.StringIO(_InPUT)
# ここからコードを記載
# 3<=n<=100
# 最小値を調べて、これが直通より小さければok

n = int(input())

c = [list(map(int, input().split())) for _ in range(n - 1)]

# c = [[0] * (n + 1) for _ in range(n + 1)]
# sum = [[0] * (n + 1) for _ in range(n + 1)]
# c = [[0] * (n)]

# for i in range(1, n):
#     values = list(map(int, input().split()))
#     for j, value in enumerate(values, start=i + 1):
#         c[i][j] = value
# print(c)
sum = [0]*(n-1)
# print(sum)
# sum[0] = c[0][0] # 駅2まで
# sum[1] = min(c[0][1],sum[0]+c[1][0]) # To3
# print(sum[1])
# sum[2] = min(sum[0]+c[1][1],sum[1]+c[2][0]) # To4
# print(sum[2])
sum[0] = c[0][0] # 駅2まで
for i in range(1,n-2+1):
    sum[i] = min(c[0][i],)
sum[1] = min(c[0][1],sum[0]+c[1][0]) # To3
print(sum[1])
sum[2] = min(c[0][2],sum[0]+c[1][1],sum[1]+c[2][0]) # To4
print(sum[2])
# sum[3] = min(c[1][4],sum[1]+c[2][3],sum[2]+c[3][2],sum[3]+c[4][1]) # To5

criterion = c[0][-1]
if criterion > sum[2]:
    print("Yes")
else:
    print("No")
# if criterion > sum[2]:
#     print("Yes")
# else:
#     print("No")
# print(c)

# itr = n-2
# for i in range(n-2):
#     for j in range(itr):
#         sum[i][j] = c[i-1][i] + c[i][j]
#     itr -= 1
# print(sum)
# max_sum = 0
# for k in range(1,n+1):
#     if k == n:
#         max_sum += c[k][-1]
#         print(c[k][-1])
#     else:
#         max_sum += c[k][k+1]
#         print(c[k][k+1])

# ans = 0
# first = list(map(int,input()))
# criterion = first[-1]
# for i in range(n-1):
#     c = list(map(int,input()))
    # if i == 0:
    #     criterion = c[-1]
    #     continue
    # else:

# all_list = [(0,0)]*n
# print(all_list)
# for i in range(n-1):
#     itr = n-1
#     for j in range(itr):
#         if itr == 1:
#             all_list[i] = int(input())
#         else:
#             all_list[i] = map(int,input().split())
#         print(all_list[i])
    
    
