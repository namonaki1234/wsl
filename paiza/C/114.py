import io
import sys
import numpy as np

# 下記に標準入力を記載
_INPUT = """\
4
apPle
error
suBway
Zb
"""
sys.stdin = io.StringIO(_INPUT)
# ここからコードを記載

N = int(input())
S = [input() for _ in range(N)]

count = 0
for i in range(N):
    
    if count == N-1 and S[N-2][len(S[N-2])-1]==S[N-1][0]:
        pass
    else:
        if S[i][len(S[i])-1]==S[i+1][0]:
            count+=1
            continue
            
        else:
            print(S[i][len(S[i])-1]+' '+S[i+1][0])
            break

    print("Yes")














