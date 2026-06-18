import io
import sys
from collections import Counter
from decimal import Decimal

# 下記に標準入力を記載
_InPUT = """\
7
3 4
6 1
5 5
2 7
7 2
1 3
4 6
"""
sys.stdin = io.StringIO(_InPUT)
# ここからコードを記載

n = int(input())

xy = []
ans = 0
for i in range(n):
    x,y = map(int,input().split())
    xy.append((x,y))

xy.sort()
# min_y = xy[0][1] + 1 # xが最小の点のy座標+1を仮の初期値としておけば1回目の処理でy_minがxが最小の点のy座標となる
min_y = float('inf') # xが最小の点のy座標+1を仮の初期値としておけば1回目の処理でy_minがxが最小の点のy座標となる
for i in range(n):
    x_i,y_i = xy[i][0],xy[i][1]
    
    if y_i < min_y:
        ans += 1
        min_y = y_i
            
print(ans)