from collections import defaultdict
import io
import sys

# 下記に標準入力を記載
_INPUT = """\
program

"""


sys.stdin = io.StringIO(_INPUT)

A= input()
B =[]
for char in A:
    if char.isupper() == True:
        B.append(char)

result = ''.join(B)
print(result)

