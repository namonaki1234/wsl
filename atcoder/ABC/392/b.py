import io
import sys
import numpy as np

# 下記に標準入力を記載
_INPUT = """\
9 1
9


"""
sys.stdin = io.StringIO(_INPUT)
# ここからコードを記載

N,M = np.array(list(map(int, input().strip().split())))
A = np.array(list(map(int, input().strip().split())))
A.sort()

count = 0
X =list(range(1,N+1))
X_result = X.copy()


for i in range(N):
    if X[i] == A[count]:
        
        X_result.remove(A[count])
        count += 1
        if count == M:
            break

        

print(len(X_result))
print(*X_result)
    
  










