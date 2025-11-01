import io
import sys


# 下記に標準入力を記載
_InPUT = """\
5
22/11
"""
sys.stdin = io.StringIO(_InPUT)
# ここからコードを記載

n = int(input())
s = input()

if len(s) % 2 == 0:
    print("No")
    exit()
else:
    mid = len(s) // 2
    if s[mid] == "/":
        for s_i in s[:mid]:
            if s_i != "1":
                print("No")
                exit()
        for s_i in s[mid+1:]:
            if s_i != "2":
                print("No")
                exit()
        print("Yes")
        exit()
    else:
        print("No")
        exit()
