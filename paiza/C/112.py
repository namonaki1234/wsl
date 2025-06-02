import io
import sys
import numpy as np

# 下記に標準入力を記載
_INPUT = """\
1
5 4 20
"""
sys.stdin = io.StringIO(_INPUT)
# ここからコードを記載

N = int(input())
s = [0] * N
f = [0] * N
t = [0] * N
for i in range(N):
    #上から順番に代入していく
    s[i], f[i], t[i] = map(int, input().split())

paiza_time = []
for i in range(N):
    paiza_time.append(s[i] + f[i]+(24-t[i]))
    
print(min(paiza_time))
print(max(paiza_time))













