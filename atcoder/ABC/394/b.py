import io
import sys
import numpy as np

# 下記に標準入力を記載
_INPUT = """\
4
cat
enate
on
c


"""
sys.stdin = io.StringIO(_INPUT)
# ここからコードを記載

N = int(input())
S = [input() for _ in range(N)]
S_value = []
for i in range(N):
    S_value.append(len(S[i]))

# for i in zip(S, S_value):
#     print(i)

for i in range(N):
    for j in range(N):
        if S_value[i] < S_value[j]:
            S[i], S[j] = S[j], S[i]
            S_value[i], S_value[j] = S_value[j], S_value[i]

print(''.join(S))
    
  










