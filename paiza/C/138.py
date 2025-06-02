import io
import sys
import numpy as np

# 下記に標準入力を記載
_INPUT = """\
5
50
50
50
60
60
"""
sys.stdin = io.StringIO(_INPUT)
# ここからコードを記載 

N = int(input())
a = [int(input()) for _ in range(N)] #sは持ってるボールの個数

count = [1 for _ in range(N)]

for i in range(N):
    for j in range(N):
        if a[i]<a[j]:
            count[i] += 1  

for i in range(N):
    print(count[i])















