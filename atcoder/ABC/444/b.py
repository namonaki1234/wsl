import io
import sys
from collections import Counter

# 下記に標準入力を記載
_InPUT = """\
99999 45
"""
sys.stdin = io.StringIO(_InPUT)
# ここからコードを記載

n,k = map(int,input().split())

n_str = str(n)
sum = 0
cnt = 0
for i in range(1,n+1):
    i = str(i)
    for j in range(len(i)):
        sum += int(i[j])
    if sum == k:
        cnt += 1
    sum = 0
        
print(cnt)


