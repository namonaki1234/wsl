import io
import sys


# 下記に標準入力を記載
_InPUT = """\
10 4
@@@.@@.@@.
"""
sys.stdin = io.StringIO(_InPUT)
# ここからコードを記載
# @のとき箱にクッキーが入っている

n,d = map(int, input().split())
s = input()

cookie_index = [i for i, c in enumerate(s) if c == "@"]

ans = n- (len(cookie_index)-d)
print(ans)
