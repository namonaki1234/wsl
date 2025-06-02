from collections import defaultdict
import io
import sys

# 下記に標準入力を記載
_INPUT = """\
400


"""


sys.stdin = io.StringIO(_INPUT)

A= input()

if A[0] =="2":
    print("Success")
else:
    print("Failure")


