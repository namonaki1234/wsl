import io
import sys

# 下記に標準入力を記載
_InPUT = """\
6 16
0 8 0 2 6 8
"""
sys.stdin = io.StringIO(_InPUT)
# ここからコードを記載

n,m = map(int,input().split())
a = list(map(int,input().split()))

a_sum = sum(a)

for i in range(n):
    a_sum_minus_a_i = a_sum - a[i]
    if a_sum_minus_a_i == m:
        print("Yes")
        exit()
print("No")
