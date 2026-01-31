import io
import sys
from collections import Counter

# 下記に標準入力を記載
_InPUT = """\
1234 12345678
"""
sys.stdin = io.StringIO(_InPUT)
# ここからコードを記載

n,k = map(int,input().split())
n_0 = n
sum = 0
for i in range(1000000000):
    sum += n
    if sum >= k:
        print(n - n_0)
        exit()
    n += 1


