from collections import defaultdict
import io
import sys

# 下記に標準入力を記載
_INPUT = """\
10
1 1 2 3 5 8 13 21 34 55


"""

sys.stdin = io.StringIO(_INPUT)

N = int(input())
A = list(map(int, input().split()))

for i in range(N):
    if i == N-1:
        break

    if A[i] < A[i+1]:
        continue
    else:
        print("No")
        exit()
print("Yes")


