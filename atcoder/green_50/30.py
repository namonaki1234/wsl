import io
import sys

# 下記に標準入力を記載
_INPUT = """\
33


"""

sys.stdin = io.StringIO(_INPUT)
# ここからコードを記載

#・全探索で解ける問題かどうか、判断できるようになる
X=int(input())

for A in range(-10**3,10**3):
    for B in range(-10**3,10**3):
        if A**5-B**5==X:
            print(A,B)
            exit()
