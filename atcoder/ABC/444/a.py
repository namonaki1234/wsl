import io
import sys

# 下記に標準入力を記載
_InPUT = """\
999
"""
sys.stdin = io.StringIO(_InPUT)
# ここからコードを記載

n = input()

cnt = 0
for i in range(3):
    if i == 2:
        break
    if n[i] == n [i+1]:
        cnt += 1

if cnt == 2:
    print("Yes")
else:
    print("No")