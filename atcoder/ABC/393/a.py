import io
import sys
import numpy as np

# 下記に標準入力を記載
_INPUT = """\
fine fine


"""
sys.stdin = io.StringIO(_INPUT)
# ここからコードを記載

A = input().strip().split()

if A[0] == 'sick' and A[1] == 'sick':
    print(1)

elif A[0] == 'sick' and A[1] == 'fine':
    print(2)

elif A[0] == 'fine' and A[1] == 'sick':
    print(3)

elif A[0] == 'fine' and A[1] == 'fine':
    print(4)
  










