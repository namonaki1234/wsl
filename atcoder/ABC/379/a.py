import io
import sys

# 下記に標準入力を記載
_InPUT = """\
379
"""
sys.stdin = io.StringIO(_InPUT)
# ここからコードを記載

n = int(input())
# aは100で割った商
a = n // 100
# bは100で割ったあまりの79を10で割った商
b = (n % 100) // 10
# cは10で割ったあまり
c = n % 10

# print(b,c,a)
print(f"{b}{c}{a} {c}{a}{b}")
# print(c,a,b)
# print(f"{c}{a}{b}")