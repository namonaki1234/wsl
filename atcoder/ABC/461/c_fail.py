import io
import sys
from collections import Counter
from decimal import Decimal

# 下記に標準入力を記載
_InPUT = """\
5 3 3
1 30
1 40
1 50
2 10
3 20
"""
sys.stdin = io.StringIO(_InPUT)
# ここからコードを記載
# n個の宝石、i番目の宝石の色はc_i(整数)、価値はv_i
# このn個の宝石の中からk個の宝石を選ぶ、ただし選んだ宝石の色がm種類以上(色で区別)
# この時選んだ宝石の価値の総和としてありうる最大値

n,k,m = map(int,input().split())
cv = []
for i in range(n):
    c,v = map(int,input().split())
    cv.append((c,v))
cv.sort(key=lambda x: x[0],reverse=True)
value_sum = 0
kind_cnt = 1
kind = cv[0][0]
skip_index = []
data = [cv[0]]
ans = cv[0][1]

j = 0
delete_index = []
visited = [False]*n
while j == n - 1:
    # ans += cv[j][0]
    # visited[j]= True
    if kind == cv[j-1][0]:
       ans = max(ans,ans+cv[j][0])
    j += 1
    
# for k in visited:
#     if visited:
#         ans += cv[j][1]

print(cv)
print(ans)
# while kind_cnt == m:
# for j in range(n):
#         if j != 0:
#             if kind == cv[j][0]:
#                 kind_cnt += 1

#         value_sum += cv[j][1]