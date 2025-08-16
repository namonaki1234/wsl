import io
import sys
from collections import defaultdict,deque,Counter
import numpy as np

# 下記に標準入力を記載
_InPUT = """\
8
greentea
"""
sys.stdin = io.StringIO(_InPUT)
# ここからコードを記載

n = int(input().strip())
s = list(input().strip())
# print(s[-3:])  # 末尾3文字を表示
if ''.join(s[-3:]) == "tea":
    print("Yes")
else:
    print("No")
