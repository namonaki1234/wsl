import io
import sys
from collections import defaultdict,deque,Counter
import numpy as np

# 下記に標準入力を記載
_InPUT = """\
2 1
2 5

"""
sys.stdin = io.StringIO(_InPUT)
# ここからコードを記載

n,k = map(int, input().split())
a = list(map(int, input().split())) 

ans = 1

for i in range(n):
    ans *= a[i]

    if ans < 10**k:
        continue
    else:
        ans = 1
    
print(ans)
    


    











