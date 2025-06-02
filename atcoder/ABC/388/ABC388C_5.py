import io
import sys
import numpy as np

# 下記に標準入力を記載
_INPUT = """\
6
2 3 4 4 7 10
"""
sys.stdin = io.StringIO(_INPUT)

N = int(input())
A = np.array(list(map(int, input().split())))

kind_count = 0

if A[N-1] / 2 < A[0]:
    print(0)
else:
    for j in range(N):
        
        half_values = A[N-1-j] / 2

        valid_indices = np.where(A <= half_values)[0]

        
        kind_count += len(valid_indices)

        
    print(kind_count)

