import io
import sys

# 下記に標準入力を記載
_InPUT = """\
1
"""
sys.stdin = io.StringIO(_InPUT)
# ここからコードを記載

n = int(input())

ans = []
for i in range(n):
    # ans.append(n)
    ans.append(str(n))
    # print(n)
    n -= 1

# print(ans)
# print(str(ans))
print(",".join(ans))