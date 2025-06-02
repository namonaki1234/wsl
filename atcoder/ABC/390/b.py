import io
import sys

# 下記に標準入力を記載
_INPUT = """\
5
3 6 12 24 48


"""
sys.stdin = io.StringIO(_INPUT)
# ここからコードを記載

N = int(input())

A = list(map(int, input().split()))

a_1 = A[0]


r = A[1]/A[0]

count = 0
a_n = a_1

for i in range(0,N-1):
    a_n *= r

    if a_n == A[i+1]:
        count += 1
    else:
        print("No")
        break
    

if count == N-1:
    print("Yes")

      