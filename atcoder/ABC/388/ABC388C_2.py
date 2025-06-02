import io
import sys

# 下記に標準入力を記載
_INPUT = """\
6
2 3 4 4 7 10

"""
sys.stdin = io.StringIO(_INPUT)
# ここからコードを記載

N = int(input())
A = list(map(int, input().split()))

kind_count = 0
if A[N-1]/2<A[0]:
    print(0)


else:
    kind_count = sum(1 for j in range(N) for i in A if A[N-1-j]/2 >= i)

    print(kind_count)




