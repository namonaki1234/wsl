import io
import sys

# 下記に標準入力を記載
_INPUT = """\
hamburg
2
cheese hamburg
"""
sys.stdin = io.StringIO(_INPUT)
# ここからコードを記載
S = input()
N = int(input())
T = input().split()


count = 0
for i in range(N):
    if S in T:
        count += 1
    
        
if count > 0:
    print("YES")
else:
    print("NO")