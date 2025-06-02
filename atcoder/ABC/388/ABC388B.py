import io
import sys

# 下記に標準入力を記載
_INPUT = """\
1 4
100 100


"""
sys.stdin = io.StringIO(_INPUT)
# ここからコードを記載

N,D = map(int,input().split())

TL = [map(int, input().split()) for _ in range(N)]
T, L = [list(i) for i in zip(*TL)]

weight = []

for i in range(D):
    for j in range(N):

    
        weight.append(T[j]*(L[j]+i+1))
    print(max(weight))



