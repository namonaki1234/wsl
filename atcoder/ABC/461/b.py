import io
import sys
from collections import Counter
from decimal import Decimal

# 下記に標準入力を記載
_InPUT = """\
5
2 4 5 1 3
4 1 5 2 3
"""
sys.stdin = io.StringIO(_InPUT)
# ここからコードを記載
# n人の木こりが一つずつ斧を持っている、全員が池に落とした
# 木こりiは自分が持っていた斧がa_iと言っている、全員本当→yes,そうでない→no
# 池の女神は斧iの持ち主がb_iと知っている

n = int(input())
a = [x for x in input().split()]
b = [x for x in input().split()]

cnt = 0
for i in range(n):
    # print(i + 1,a[int(b[i]) - 1])
    if (i + 1) == int(a[int(b[i]) - 1]):
        cnt += 1

if cnt == n:
    print("Yes")
else:
    print("No")
