import io
import sys
from collections import defaultdict,deque,Counter
import numpy as np

# 下記に標準入力を記載
_InPUT = """\
atcoder


"""
sys.stdin = io.StringIO(_InPUT)
# ここからコードを記載

s = list(input().strip())

english_words = ["red", "blue", "green"]
atcoder_words = ["SSS", "FFF", "MMM"]

# print(''.join(s))
if ''.join(s) in english_words:
    if ''.join(s) == "red":
        print("SSS")
    elif ''.join(s) == "blue":
        print("FFF")
    elif ''.join(s) == "green":
        print("MMM")
else:
    print("Unknown")