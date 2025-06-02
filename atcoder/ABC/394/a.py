from collections import defaultdict
import io
import sys

# 下記に標準入力を記載
_INPUT = """\
22222222222
22222222222
"""

sys.stdin = io.StringIO(_INPUT)

S = input()
S_list = []

for char in S:
    S_list.append(char)

length = list(range(len(S)))
length.sort(reverse = True)

for i in length:
    if S[i] != '2':
        S_list.pop(i)
        

ans = ''.join(S_list)
print(ans)


