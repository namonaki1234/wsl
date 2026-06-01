import io
import sys

# 下記に標準入力を記載
_InPUT = """\
5
"""
sys.stdin = io.StringIO(_InPUT)
# ここからコードを記載

x = int(input())
hw = list("HelloWorld")
# print(hw)
hw.pop(x-1)


print("".join(hw))