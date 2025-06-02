from collections import defaultdict
import io
import sys
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')
sys.stderr = io.TextIOWrapper(sys.stderr.buffer, encoding='utf-8')

# 下記に標準入力を記載
_INPUT = """\
8
chokudai
chokudai

"""


sys.stdin = io.StringIO(_INPUT)

N= int(input())
S = input().rstrip()
T = input().rstrip()

ans = 0
for i in range(N):
    if S[i] != T[i]:
        ans += 1
        
print(ans)




