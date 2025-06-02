from collections import defaultdict
import io
import sys
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')
sys.stderr = io.TextIOWrapper(sys.stderr.buffer, encoding='utf-8')

# 下記に標準入力を記載
_INPUT = """\
4
"""


sys.stdin = io.StringIO(_INPUT)

N= int(input())

S =["-"]*N
for i in range(N):
    # 偶数の時
    if N % 2 == 0:
        if i == N/2-1:
            S[i] = "="
            S[i+1] = "="
    # 奇数の時
    else:
        if i == N//2:
            S[i] = "="
   
print("".join(S))
# print(7//2)
  




