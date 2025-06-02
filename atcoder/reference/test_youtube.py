import io
import sys
from itertools import product,combinations
from heapq import heapify,heappop,heappush

# 辞書、リスト、タプルを進化させたものが入っているライブラリ
from collections import deque,defaultdict
# 下記に標準入力を記載
_INPUT = """\

"""
sys.stdin = io.StringIO(_INPUT)
# ここからコードを記載

li_0 = list(combinations(range(1,5),2))
# print(li_0)

# 直積はbit全探索で使う
N = 3
li_1 = product([0,1],repeat=N)

# print(list(li_1))

# dequeは先頭と末尾に対しての追加や削除が高速
li_2 = [0,1,2,3,4,5]
q_0 = deque(li_2)
# print(q)

# while q:
#     l = q.popleft()
#     r = q.pop()
#     print(l,r,q)


# defaultdictは辞書の初期化を省略できる
d = defaultdict(int)
d[0] = 1
d[1] = 2
d[2] = 3
# print(d)

# ヒープキュー（優先度付きキュー）は最小値,最大値を取り出したり、追加したりするのが高速
# 何回も要素を追加するたびに最小値を取り出すときに使う
q_1 = [3,1,4,1,0,6,5,3,5]
# heapify(q_1)
# print(q_1)

# while q_1:
#     x = heappop(q_1)
#     print(x,q_1)

# if文の分岐後のprint("Yes")の後のexit()で余計なprintをしないようにする

# for文はrange関数を使うときに、range(0,len(q_1),2)のようにすると、偶数番目の要素だけを取り出すことができる
# for i in range(0,len(q_1),2):
#     print(q_1[i])

# 最大値を求める問題では最初に最小値でansを定義して、それより大きければ更新するという方法で計算量N回で実行できる

tst = [-1,3,4,-6]

ans = 0
for i in tst:
    # print(abs(i))   
    ans = max(ans,abs(i))
print(ans)



