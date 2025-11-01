import io
import sys


# 下記に標準入力を記載
_InPUT = """\
10 b a
acaabcabba
"""
sys.stdin = io.StringIO(_InPUT)
# ここからコードを記載

n, c1, c2 = input().split()
s = input()

target_index = [i for i, c in enumerate(s) if c == c1]
# print(target_index)
ans = [c1 if i in target_index else c2 for i in range(int(n))]
print("".join(ans))