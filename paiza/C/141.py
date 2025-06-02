import io
import sys
import numpy as np

# 下記に標準入力を記載
_INPUT = """\
5
kato
kato
sato
sato
sato
"""
sys.stdin = io.StringIO(_INPUT)
# ここからコードを記載

N = int(input())
a = [input() for _ in range(N)] #sは持ってるボールの個数

print(max(a, key=a.count))















