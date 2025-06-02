from collections import defaultdict
import io
import sys

# 下記に標準入力を記載
_INPUT = """\
13
1 1 1 1 1 1 1 1 1 1 1 1 1


"""

sys.stdin = io.StringIO(_INPUT)

N = int(input())
A = list(map(int, input().split()))

for i in range(N):
    if i == N-1:
        break
    if A[i] == A[i+1] and A[i] == A[i+2]:
        print("Yes")
        exit()

print("No")


