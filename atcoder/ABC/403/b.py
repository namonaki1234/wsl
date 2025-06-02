import io
import sys
from collections import defaultdict,deque,Counter

# 下記に標準入力を記載
_INPUT = """\
tak??a?h?
nashi

"""
sys.stdin = io.StringIO(_INPUT)
# ここからコードを記載

T = input()
U = input()
indices = [i for i in range(len(T)) if T[i] == "?"]
alphabets = [chr(ord("a") + i) for i in range(26)]
for a in alphabets:
    for b in alphabets:
        for c in alphabets:
            for d in alphabets:
                S = list(T)
                S[indices[0]] = a
                S[indices[1]] = b
                S[indices[2]] = c
                S[indices[3]] = d
                S = "".join(S)
                if U in S:
                    print("Yes")
                    exit()
print("No")
     





    


    











