import io
import sys
import numpy as np

# 下記に標準入力を記載
_INPUT = """\
6
2 3 4 4 7 10

"""
sys.stdin = io.StringIO(_INPUT)
# ここからコードを記載

N = int(input())
A = np.array(list(map(int, input().split())))

kind_count = 0
if A[N-1]/2<A[0]:
    print(0)


else:
    for j in range(N):
        if np.sum(A <= A[N-1-j]/2) == 0:  
            break  
        kind_count += np.sum(A <= A[N-1-j]/2)


    print(kind_count)





