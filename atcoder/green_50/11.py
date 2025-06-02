import io
import sys

# 下記に標準入力を記載
_INPUT = """\
2 3
2 1
5 10

"""

sys.stdin = io.StringIO(_INPUT)
# ここからコードを記載

# 目的はTLE回避
N,K=map(int, input().split())

friends=[]
for i in range(N):
    A,B=map(int, input().split())
    friends.append([A,B])

friends.sort()

now_village=0

now_village+=K

for i in range(N):
    friend_village=friends[i][0]
    friend_money=friends[i][1]

    if friend_village<=now_village:
        now_village+=friend_money
    else:
        break

print(now_village)
