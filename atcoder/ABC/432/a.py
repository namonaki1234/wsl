import io
import sys


# 下記に標準入力を記載
_InPUT = """\
9 1 9
"""
sys.stdin = io.StringIO(_InPUT)
# ここからコードを記載
a, b, c = map(int, input().split())

list = [a, b, c]
list.sort(reverse=True)
# print(list)

print(''.join(map(str, list)))
