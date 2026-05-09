import io
import sys

# 下記に標準入力を記載
_InPUT = """\
5
1 2 3 4 5
3
"""
sys.stdin = io.StringIO(_InPUT)
# ここからコードを記載

n = int(input())
a = [x for x in input().split()]
x = int(input())

print(a[x-1])