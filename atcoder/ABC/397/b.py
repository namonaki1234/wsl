import io
import sys
from collections import defaultdict,deque

# 下記に標準入力を記載
_INPUT = """\
oiii


"""
sys.stdin = io.StringIO(_INPUT)
# ここからコードを記載

S = input()

ans = 0
if S[0] == "o":    
        ans += 1


for j in range(len(S)):
    
    if S[j] == "i" and j+1 < len(S) and S[j+1] == "i":
        ans += 1
    elif S[j] == "i" and j+1 < len(S) and S[j+1] == "o":
        continue
    elif S[j] == "o" and j+1 < len(S) and S[j+1] == "o":
        ans += 1
    elif S[j] == "o" and j+1 < len(S) and S[j+1] == "i":
        continue
    elif S[j] == "i" and j== len(S)-1:
        ans += 1
    
print(ans)











    











