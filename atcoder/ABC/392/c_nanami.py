import io
import sys
import numpy as np
from collections import deque
import bisect



# 下記に標準入力を記載
_INPUT = """\
4
4 3 2 1
2 3 1 4


"""
sys.stdin = io.StringIO(_INPUT)
# ここからコードを記載
# N,M = np.array(list(map(int, input().strip().split())))
# P[i]が見つめている人、P[i]がつけているゼッケンがQ[i]、i=1からNまでの順番に並んでいる

N = int(input())
P_view = np.array(deque(map(int, input().strip().split()))) 
Q = np.array(deque(map(int, input().strip().split())))


P = deque(range(1,N+1))

result = []

# bisect.bisect_left(Q, i) ==
QP = []
for i in range(1,N+1):
    QP.append((Q[i-1],P_view[i-1]))
    QP_sorted = sorted(QP)
   
    
for i in range(N):
    P_index = QP_sorted[i][1]-1
    result.append(QP[P_index][0])




print(*result)


















