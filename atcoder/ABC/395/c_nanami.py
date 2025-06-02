import io
import sys
import re

# 下記に標準入力を記載
_INPUT = """\
WACWA


"""
sys.stdin = io.StringIO(_INPUT)
# ここからコードを記載

S = input()
length = list(range(len(S)))
length.sort(reverse=True)

while re.search(r'W+A', S) != None:
    replace = lambda m: 'A' + 'C' * (len(m.group(0)) - 1)
    # for i in length:
        # if i+1 >= len(S):
        #     break
        
    # if r'W+A' in S:
    S = re.sub(r'W+A', replace, S)
        # match S[:i]:
        #     case r'W+A':
                # S = re.sub(r'W+A', replace, S)
                # S[:i] + replace + S[i+len(replace):]
            # case 'WA':
            #     S = S[:i] + 'AC' + S[i+2:]
            # case _:
            #     pass


print(S)
    
  










