import io
import sys
from collections import Counter


# 下記に標準入力を記載
_InPUT = """\
500000
"""
sys.stdin = io.StringIO(_InPUT)
# ここからコードを記載

n = int(input())

counter_n = Counter(str(n))

if counter_n["1"] == 1 and counter_n["2"] == 2 and counter_n["3"] == 3:
    print("Yes")
else:
    print("No")