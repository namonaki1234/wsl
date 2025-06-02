import io
import sys

# 下記に標準入力を記載
_INPUT = """\
3 3 10
60 2 2 4
70 8 7 9
50 2 3 9


"""

sys.stdin = io.StringIO(_INPUT)
# ここからコードを記載

# Bit全探索を実装できるようになる
N,M,X=map(int, input().split())

C = []  # 各本の値段
A = []  # 各本のアルゴリズムへの効果, i番目の本を読むと j番目のアルゴリズムがどれだけ上がるか
for i in range(N):
    CA = list(map(int, input().split()))
    C.append(CA[0])
    A.append(CA[1:])


ans=10**20
for i in range(1<<N):
    cost=0
    skill=[0]*M

    for shift in range(N):
        if i>>shift & 1==1:
            cost+=C[shift]
            for j in range(M):
                skill[j]+=A[shift][j]

    if X<=min(skill):
        ans=min(ans,cost)

if ans==10**20:
    print(-1)
else:
    print(ans)
