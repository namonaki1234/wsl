import io
import sys

# 下記に標準入力を記載
_InPUT = """\
4
4 3 2 5
"""
sys.stdin = io.StringIO(_InPUT)
# ここからコードを記載

n = int(input())
a = list(map(int, input().split()))
for i in range(n):
    ans = -1
    for j in range(i - 1, -1, -1):
        if a[j] > a[i]:
            ans = j + 1
            break
    print(ans)
