import io
import sys
import numpy as np

# 下記に標準入力を記載
_INPUT = """\
5
0
20
10
15
2
3
1 3 5
3 2 3
2 1 10
"""
sys.stdin = io.StringIO(_INPUT)
# ここからコードを記載 

N = int(input())
s = [int(input()) for _ in range(N)] #sは持ってるボールの個数
M = int(input())
a = [0] * M
b = [0] * M
x = [0] * M
for i in range(M):
    #上から順番に代入していく
    a[i], b[i], x[i] = map(int, input().split())

for i in range(M):
    if x[i] <= s[a[i]-1]:
        s[a[i]-1]-=x[i]
        s[b[i]-1]+=x[i]
    else:
        s[b[i]-1]+=s[a[i]-1]
        s[a[i]-1]-=s[a[i]-1]
        

for i in range(N):
    print(s[i])















