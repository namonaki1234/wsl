import io
import sys

# 下記に標準入力を記載
_INPUT = """\
2
FLIP
2
2 0 0
1 1 4


"""

sys.stdin = io.StringIO(_INPUT)
# ここからコードを記載

# 目的は反転、並び替えの操作における状態管理変数による計算量削減
# 実際に入れ替え操作をする代わりに、今『通常状態』なのか『入替状態』なのかを変数で管理する、という方法を使う。
N=int(input())
S=input()
Q=int(input())

S="0"+S
S=list(S)

flip=False

for i in range(Q):
    T,A,B=map(int, input().split())

    if T==1:
        if flip==False:
            S[A],S[B]=S[B],S[A]

        else:
            if A<=N:
                A+=N
            else:
                A-=N
            if B<=N:
                B+=N
            else:
                B-=N
            S[A],S[B]=S[B],S[A]

    # T=2
    else:
        if flip==True:
            flip=False
        else:
            flip=True

S_left=S[1:N+1]
S_right=S[N+1:]

if flip==False:
    ans=S_left+S_right
else:
    ans=S_right+S_left

print("".join(ans))
