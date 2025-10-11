import io
import sys


# 下記に標準入力を記載
_InPUT = """\
ATCODER

"""
sys.stdin = io.StringIO(_InPUT)
# ここからコードを記載

s = input()

mid_index = (len(s) + 1) // 2 - 1

s_list = list(s)
s_list[mid_index] = ""
print("".join(s_list))