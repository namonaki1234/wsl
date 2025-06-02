import io
import sys
import numpy as np

# 下記に標準入力を記載
_INPUT = """\
3 3 9



"""
sys.stdin = io.StringIO(_INPUT)
# ここからコードを記載

A = np.array(list(map(int, input().strip().split())))

if A[0]*A[1]==A[2] or  A[0]*A[2]==A[1] \
or A[1]*A[0]==A[2] or A[1]*A[2]==A[0] \
or A[2]*A[0]==A[1] or A[2]*A[1]==A[0] :
    print("Yes")
else:
    print("No")
  










