import io
import sys

# 下記に標準入力を記載
_INPUT = """\
3 4
2
1

"""
sys.stdin = io.StringIO(_INPUT)
# ここからコードを記載
N,D = map(int,input().split())
d = [int(input().strip()) for _ in range(N-1)]
x = D 

for i in range(N-1):
    x += int(D-d[i])

print(D*x)
