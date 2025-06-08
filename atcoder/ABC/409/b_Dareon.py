import io
import sys
from collections import defaultdict,deque,Counter
import numpy as np

# 下記に標準入力を記載
_InPUT = """\
7
1 6 2 10 2 3 2
"""
sys.stdin = io.StringIO(_InPUT)
# ここからコードを記載

N=int(input())
A=list(map(int,input().split()))
A.sort()
ans=0
L=[0]*200
L[0]=N
for i in range(N):
    if A[i]>150:
        continue
    L[A[i]+1]-=1
for i in range(1,200):
    L[i]+=L[i-1]
for i in range(200):
    if L[i]>=i:
        ans=max(ans,i)
print(ans)
