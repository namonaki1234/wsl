import io
import sys

# 下記に標準入力を記載
_INPUT = """\
5 2
8 7 6
rsrpr

"""

sys.stdin = io.StringIO(_INPUT)
# ここからコードを記載

#・貪欲を理解し、実装できるようになる
N,K=map(int, input().split())
R,S,P=map(int, input().split())
T=input()

hands=[]
ans=0

for i in range(N):
    if T[i]=="r":
        if i<K:
            ans+=P
            hands.append("p")
        elif K<=i and hands[i-K]!="p":
            ans+=P
            hands.append("p")
        else:
            hands.append("x")

    if T[i]=="s":
        if i<K:
            ans+=R
            hands.append("r")           
        elif K<=i and hands[i-K]!="r":
            ans+=R
            hands.append("r")
        else:
            hands.append("x")

    if T[i]=="p":
        if i<K:
            ans+=S
            hands.append("s")          
        elif K<=i and hands[i-K]!="s":
            ans+=S
            hands.append("s")
        else:
            hands.append("x")

print(ans)
