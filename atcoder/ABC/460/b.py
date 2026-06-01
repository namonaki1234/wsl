import io
import sys
from collections import Counter
from decimal import Decimal

# 下記に標準入力を記載
_InPUT = """\
7
0 0 2 2 3 2
0 0 2 2 3 1
1 2 5 3 2 1
5 4 2 8 8 3
2 1 5 5 1 2
0 0 1 0 0 1
0 0 500000000 1 1000000000 500000000
"""
sys.stdin = io.StringIO(_InPUT)
# ここからコードを記載

t = int(input())

for i in range(t):
    x1,y1,r1,x2,y2,r2 = map(int,input().split())
    a = Decimal(str((x2 - x1)**2))
    b = Decimal(str(((y2 - y1)**2)))
    d = (( a + b ) ** Decimal(str(0.5)))
    # print(d)
    R_plus = r1 + r2
    # print(R_plus)
    R_minus = abs(r2 - r1)
    # print(R_minus)
    if R_minus <= d <= R_plus:
        print("Yes")
    else:
        print("No")

