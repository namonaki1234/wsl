import io
import sys
import numpy as np

# 下記に標準入力を記載
_INPUT = """\
--2000--01---01--Fizz----Buzz--
"""
sys.stdin = io.StringIO(_INPUT)
# ここからコードを記載

S = list(input())


count = 0
minus_index = []
for i in range(len(S)):
    if S[i] == "-":
        if i >= len(S)-1 and S[len(S)-1] == "-":
            pass
        else:
            if S[i+1] == "-":
                minus_index.append(i+1)
            
            

for i in sorted(minus_index,reverse=True):
    
    S.pop(i)            
    

print("".join(S))














