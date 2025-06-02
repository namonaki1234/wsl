import io
import sys
import numpy as np

# 下記に標準入力を記載
_INPUT = """\
ALANTURING
"""
sys.stdin = io.StringIO(_INPUT)
# ここからコードを記載

S = input()

before_s = ["A","E","G","I","O","S","Z"]
after_s = [4,3,6,1,0,5,2]

count = 0
for i in range(len(S)):
    for j in range(len(before_s)):
        if  before_s[j]  in S:
            S = S.replace(before_s[j],str(after_s[j]))


print(S)














