import io
import sys
from collections import Counter

# 下記に標準入力を記載
_InPUT = """\
2
algorithm heuristic
"""
sys.stdin = io.StringIO(_InPUT)
# ここからコードを記載

n = int(input())
s = [x for x in input().split()]
c= []
for i in range(n):
    if s[i][0] in ["a","b","c"]:
        c.append("2")
    elif s[i][0] in ["d","e","f"]:
        c.append("3")
    elif s[i][0] in ["g","h","i"]:
        c.append("4")
    elif s[i][0] in ["j","k","l"]:
        c.append("5")
    elif s[i][0] in ["m","n","o"]:
        c.append("6")
    elif s[i][0] in ["p","q","r","s"]:
        c.append("7")
    elif s[i][0] in ["t","u","v"]:
        c.append("8")
    elif s[i][0] in ["w","x","y","z"]:
        c.append("9")

print("".join(c))
# print("".join(list(str(c))))
# print((str(c)))

