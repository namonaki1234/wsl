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

A =[1]*(K)

if N >= K:
    q_A = deque(A)
    for i in range(K+1):
        plus = sum(q_A)
        q_A.popleft()
        q_A.append(plus)

    print(q_A[-1]%10**9)
else:
    print(1)






