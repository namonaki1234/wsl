import io
import sys

# 下記に標準入力を記載
_INPUT = """\
4 2 3 1 5

"""
sys.stdin = io.StringIO(_INPUT)
# ここからコードを記載

A = list(map(int, input().split()))

sorted_A = sorted(A)


change_number =[]

for i in range(len(A)):
    if A[i] != sorted_A[i]:
        change_number.append(A[i])

if len(change_number) == 2 and abs(change_number[0] - change_number[1]) == 1:
    print("Yes")
else:
    print("No")
        


