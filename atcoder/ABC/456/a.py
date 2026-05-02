import io
import sys

# 下記に標準入力を記載
_InPUT = """\
15
"""
sys.stdin = io.StringIO(_InPUT)
# ここからコードを記載
# サイコロ3つ、出ためがxになるか判定

x = int(input())

for i in range(1,7):
    for j in range(1,7):
        for k in range(1,7):
            if (i + j + k) == x:
                print("Yes")
                exit()

print("No")