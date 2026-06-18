import io
import sys
from collections import Counter
from decimal import Decimal

# 下記に標準入力を記載
_InPUT = """\
4
1 2
1 3
1 2
3 1 2 3
"""
sys.stdin = io.StringIO(_InPUT)
# ここからコードを記載
# 人1~nが互いにギフトを送り合う
# 人iは人a_i,1・・人a_i,k_iのk_i人にギフトを送った
# i=1,..,nに対し、人iにギフトを送った人を全て求める
# 出力：x b_1 .. b_x
# xは人iにギフトを送った人数,bはその人の番号


n = int(input())
# k_a = [0]*n 

# for i in range(n):
#     k_a[i] = [int(x) for x in input().split()]

# print(k_a)
# ans = [[0]]*(n+1)
ans = [[] for _ in range(n+1)]

for i in range(1,n+1):
    k_a = [int(x) for x in input().split()]
    k = k_a[0]
    a = k_a[1:]
    for a_j in a:
        ans[a_j].append(i)

# print(ans)
for j in range(len(ans)):
    if j == 0:
        continue
    print(len(ans[j]),*ans[j])