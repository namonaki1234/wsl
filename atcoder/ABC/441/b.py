import io
import sys
from collections import Counter

# 下記に標準入力を記載
_InPUT = """\
13 11
defghiqsvwxyz
acejmoqrtwx
15
qhsqzhd
jcareec
wwqxqew
wxqxwex
jxxrtwa
trtqjxe
sqyggse
xxqwxew
xewwxxw
wwqxwex
xqqxqwq
qxxexxe
teqeroc
eeeqqee
vxdevyy
"""
sys.stdin = io.StringIO(_InPUT)
# ここからコードを記載

n,m = map(int,input().split())
# s = list(input())
s = input()
t = input()
# print(s,t)
# t = list(input())
# print(s,t)
q = int(input())
# words = set()
for i in range(q):
    w = set(input())
    w = list(w)
    w.sort()
    # w = str(w)
    w = "".join(w)
    it_s = iter(s)
    it_t = iter(t)
    # print(w)
    if all(ch in it_s for ch in w) and all(ch in it_t for ch in w):
        print("Unknown")
        continue
    it_s = iter(s)
    it_t = iter(t)
    if all(ch in it_s for ch in w):
        print("Takahashi")
    elif all(ch in it_t for ch in w):
        print("Aoki")
    else:
        print("Unknown")

    # if w in s:
    #     print("Takahashi")
    # elif w in t:
    #     print("Aoki")
    # else:
    #     print("Unknown")

# print(words)
# for word in words:
    
# unique_words =set()
# for word in words:
#     unique_words.add(word)

# print(unique_words)
