import io
import sys

# 下記に標準入力を記載
_INPUT = """\
2 2000
450 210
"""
sys.stdin = io.StringIO(_INPUT)
# ここからコードを記載
N,L = map(int,input().split())
T = list(map(int,input().split()))

sum = 0
if max(T) < L:
    for i in range(N):
        sum += T[i]
    print(sum)
else:
    for i in range(N-1):
        T_sorted = sorted(T)
        sum += T_sorted[i]
    
    half_max = max(T)/2

    sum += half_max
    print(int(sum))
        

