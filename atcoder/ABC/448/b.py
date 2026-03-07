import io
import sys
from collections import Counter

# 下記に標準入力を記載
_InPUT = """\
15 10
7 94 100 82 63 81 75 2 76 73
10 44
5 77
10 47
7 32
2 82
5 90
3 37
6 70
6 28
3 25
2 26
10 56
1 69
5 46
7 26
"""
sys.stdin = io.StringIO(_InPUT)
# ここからコードを記載
# m種類のコショウ、それぞれc_i個ある
# 高橋君はn個の注文をしたが、相性の都合でi番目の料理には種類a_iのコショウしか使えない(上限b_iグラム)
# 料理にかけたコショウのグラムを最大にしたい、その最大を出力

n,m = map(int,input().split())
c = list(map(int,input().split()))

# a = [0]*n
# b = [0]*n
# for i in range(n):
#     a[i],b[i] = map(int,input().split())
# print(a,b)

sum = 0
for i in range(n):
    a,b = map(int,input().split())
    if b < c[a-1]:
        sum += b
        c[a-1] -= b
    else:
         sum += c[a-1]
         c[a-1] = 0

print(sum)