import io
import sys

# 下記に標準入力を記載
_InPUT = """\
48 25
"""
sys.stdin = io.StringIO(_InPUT)
# ここからコードを記載
# n個の座席、m人、一つ飛ばしに座れるか判定

n,m = map(int,input().split())

if n%2 == 0 and n>=2*m:
    print("Yes")
    exit()
elif n%2 != 0:
    n_half = (n+1)/2
    if m<= n_half:
        print("Yes")
        exit()

print("No")
