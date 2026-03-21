import io
import sys

# 下記に標準入力を記載
_InPUT = """\
4
25 40 65
30 55
25
"""
# 5
# 25 40 65 80
# 30 55 60
# 25 30
# 20
# """
sys.stdin = io.StringIO(_InPUT)
# ここからコードを記載
# 3<=n<=100

n=int(input())
C=[[None]*(i+1)+list(map(int,input().split())) for i in range(n-1)]
# print(C)
ans="No"
for a in range(n):
    for b in range(a+1,n):
        for c in range(b+1,n):
            if C[a][b]+C[b][c]<C[a][c]:
                ans="Yes"
print(ans)
