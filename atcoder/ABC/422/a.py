import io
import sys
from collections import defaultdict,deque,Counter

# 下記に標準入力を記載
_InPUT = """\
3-3

"""
sys.stdin = io.StringIO(_InPUT)
# ここからコードを記載

s = input()

if s[-1] == "8":
    s_list = list(s)
    s_list[-1] = "1"
    s_list[0] = str(int(s_list[0]) + 1)
    s = "".join(s_list)
else:
    s_list = list(s)
    s_list[-1] = str(int(s_list[-1]) + 1)
    s = "".join(s_list)
print(s)