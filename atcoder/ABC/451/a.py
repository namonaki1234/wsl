import io
import sys

# 下記に標準入力を記載
_InPUT = """\
illegal
"""
sys.stdin = io.StringIO(_InPUT)
# ここからコードを記載

s = input()
if len(s)%5==0:
    print("Yes")
else:
    print("No")    

