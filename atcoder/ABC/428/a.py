import io
import sys


# 下記に標準入力を記載
_InPUT = """\
1 1 666 428

"""
sys.stdin = io.StringIO(_InPUT)
# ここからコードを記載

s,a,b,x = map(int, input().split())

ans = 0
while x > 0:
    if x <= a:
        ans += s*x
        break
    x -= a
    ans += s*a
    x -= b

print(ans)