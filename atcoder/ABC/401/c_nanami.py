import io
import sys
from collections import defaultdict,deque,Counter
import math
# from more_itertools import chunked

# 下記に標準入力を記載
_INPUT = """\
1000000 500000
"""
sys.stdin = io.StringIO(_INPUT)
# ここからコードを記載
N, K = map(int, input().split())

A =[1]*(K+1)

if N >= K:
    q_A = deque(A)
    for i in range(K,N):
        if K<=i<2*K:
            A[K] = K
            A[i+1] = A[2*i-1]
        elif i >= 2*K:
            A[i] = 2*A[i-1]-A[(i-1)/2]
            
    print(A[N]%10**9)
else:
    print(1)



