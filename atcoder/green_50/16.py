import io
import sys

# 下記に標準入力を記載
_INPUT = """\
4 330
0 1 10 100
1 0 20 200
10 20 0 300
100 200 300 0

"""

sys.stdin = io.StringIO(_INPUT)
# ここからコードを記載

# 目的はitertoolsを使えるようになる
N,K=map(int, input().split())

time=[]

for i in range(N):
    T=list(map(int, input().split()))
    time.append(T)

ans=0

import itertools
for root in itertools.permutations(range(1,N)):
    now_time=0
    now_time+=time[0][root[0]]
    now_city=root[0]

    for i in range(1,N-1):
        to_city=root[i]
        now_time+=time[now_city][to_city]
        now_city=to_city

    now_time+=time[now_city][0]
    if now_time==K:
        ans+=1

print(ans)


"""
それぞれの書き方を以下に記載するのでコンテスト中にコピペできるよう保存しておこう。
import itertools
# 順列
#(1,2,3),(1,3,2),(2,1,3),(2,3,1),...,(3,2,1)
for seq in itertools.permutations(range(1,4)):
# 重複なしの組み合わせ
# (1,2,3),(1,2,4),...,(7,8,9)
for seq in itertools.combinations(range(1,10),3):
# 重複ありの組み合わせ
#(1,1,1),(1,1,2),...,(9,9,9)
for seq in itertools.combinations_with_replacement(range(1,10),3):
# 直積
#(1,1),(1,2),(1,3),(2,1),(2,2),...,(3,3)
for seq in itertools.product(range(1,4),range(1,4)):
"""