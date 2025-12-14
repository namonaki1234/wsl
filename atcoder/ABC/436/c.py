import io
import sys
from collections import Counter

# 下記に標準入力を記載
_InPUT = """\
4 3
1 1
2 2
2 3
"""
sys.stdin = io.StringIO(_InPUT)
# ここからコードを記載

n,m = map(int,input().split())

R,C = [],[]
for i in range(m):
    r,c = map(int,input().split())
    R.append(r)
    C.append(c)


print(R,C)
