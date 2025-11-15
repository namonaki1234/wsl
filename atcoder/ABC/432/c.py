import io
import sys

# 下記に標準入力を記載
_InPUT = """\
3 6 8
11 10 13
"""
sys.stdin = io.StringIO(_InPUT)
# ここからコードを記載
# 小さい飴をがxグラム、大きい飴がyグラム、ただしx < y
# n人の子供に飴を配るとき、各子供がもらう飴の重さの合計がすべて同じになるように配ることができるか？

n,x,y = map(int, input().split())
a = list(map(int, input().split()))

# weight_limit = y * min(a)

sum_pair = [[] for _ in range(n)]
# sum_pair = [[]]*n
for i in range(n):
    for j in range(a[i]+1):
        sum_pair[i].append((j*x + (a[i]-j)*y, j, a[i]-j))
        # sum_pair[i].append(j*x + (a[i]-j)*y, j, a[i]-j)

print(sum_pair)

# for sum_list in sum_pair:
#     for sum_val in sum_list:
        # print(sum_val)


























# dp =[max(a)][max(a)*y]
# weight_limit =  y * max(a)
# small_count = 0
# large_count = 0
# dp = [[False]*(y*max(a)+1) for i in range(n)]
# # print(dp)

# dp[0][0] = True

# for i in range(n):
#     for j in range(y*a[i]+1):
#         if j < a[i]:
#             if dp[i-1][j]:
#                 dp[i][j] = True
#         else:
#             if dp[i-1][j] or dp[i-1][j - a[i]] :
#                 dp[i][j] = True
#             else:
#                 dp[i][j] = False

# print(dp)

# if dp[n][y*max(a)]:
#     print("Yes")
# else:
#     print("No")

        

