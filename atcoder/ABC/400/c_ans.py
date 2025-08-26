import io
import sys
from math import isqrt
# from more_itertools import chunked

# 下記に標準入力を記載
_INPUT = """\
20

"""
sys.stdin = io.StringIO(_INPUT)
# ここからコードを記載

"""
正の整数(a,b) を用いて2^a x b^2  と表せる整数は、a が偶数なら平方数を4 倍した数、a が奇数なら平方数を2 倍した数になります。
よって、N 以下の良い整数の個数は⌊√N/4⌋+⌊√N/2⌋ 個になります。
実装の際に sqrt 関数をそのまま使うと小数点誤差により正しい値が得られない場合があるので注意してください。
"""

n = int(input())
print(isqrt(n // 2) + isqrt(n // 4))
