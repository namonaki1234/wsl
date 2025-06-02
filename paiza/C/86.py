import io
import sys
import numpy as np

# 下記に標準入力を記載
_INPUT = """\
Torvalds

"""
sys.stdin = io.StringIO(_INPUT)
# ここからコードを記載

S = list(input())



count = 0

A = ["a", "e", "i", "o", "u", "A", "E", "I", "O", "U"]

S_result = S.copy()
for i in range(len(S)):
    for j in range(len(A)):
        if S[i] == A[j]:
            
            S_result.remove(A[j])
            
        
        
S_result = ''.join(S_result)
print(S_result  )
    
  










