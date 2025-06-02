import io
import sys

# 下記に標準入力を記載
_INPUT = """\
WWA


"""
sys.stdin = io.StringIO(_INPUT)
# ここからコードを記載

S = input()

def match(S):
    for i in range(len(S)):
        if i+1 >= len(S):
            break
        double_string = S[i] + S[i+1]
        
        match double_string:
            case 'WA':
                S = S[:i] + 'AC' + S[i+2:]
            case _:
                pass
    return S

while 'WA' in S:
    S = match(S)

print(S)
    
  










