import io
import sys
from collections import defaultdict,deque

# 下記に標準入力を記載
_INPUT = """\
5
2
2
2
2
2


"""
sys.stdin = io.StringIO(_INPUT)
# ここからコードを記載

Q = int(input())
q = [tuple(map(int,input().split())) for _ in range(Q)]
A = [0] * (10**2)

A_q = deque(A)

for q_i in range(Q):
    if q[q_i][0] == 1:
        A_q.append(q[q_i][1])
    elif q[q_i][0] == 2:
        print(A_q.pop())













    











