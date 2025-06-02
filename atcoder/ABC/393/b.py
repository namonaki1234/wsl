import io
import sys
import numpy as np

# 下記に標準入力を記載
_INPUT = """\
AABAAABBAEDCCCD

"""
sys.stdin = io.StringIO(_INPUT)
# ここからコードを記載

S = input().strip()
S = '0'+S
A_index = []
B_index = []
C_index = []

for i in range(len(S)):
    if S[i] == 'A':
        A_index.append(i)
    if S[i] == 'B':
        B_index.append(i)
    if S[i] == 'C':
        C_index.append(i)

count = 0
for i in range(len(A_index)):
    for j in range(len(B_index)):
        for k in range(len(C_index)):
            if A_index[i] < B_index[j] < C_index[k] and (B_index[j] - A_index[i]) == (C_index[k] - B_index[j]):
                count += 1



print(count)
        


  










