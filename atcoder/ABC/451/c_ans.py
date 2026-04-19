import io
import sys
from collections import Counter

# 下記に標準入力を記載
_InPUT = """\
5
1 5
1 7
1 8
2 7
1 3
"""
sys.stdin = io.StringIO(_InPUT)
# ここからコードを記載

# q = int(input())
# query = [list(map(int, input().split())) for _ in range(q)]

#1は庭に高さh の木を新たに1 本追加する
#2はいま庭にある木のうち、高さがh 以下の木をすべて削除する

"""
2. ふつうにやると何が困るのか

たとえば木をただのリストで持って、
trees = [5, 7, 8]とする。

クエリ2で7が来たら「7以下を全部削除」したいが、そのためには
5 は削除、7 は削除、8 は残す
と、全部見ないといけない

クエリが最大 3×10^5 個あるので、毎回全部見ると遅くなる

3. この問題のコツ

この問題では「h 以下の木を全部削除する」のが目的
つまり大事なのは、いちばん小さい木が何か

なぜなら、
いちばん小さい木が h 以下なら削除する
その次に小さい木も h 以下なら削除する
そうやって、もう小さい木が h 以下でなくなるまで続ければよい

つまり、常に最小値をすぐ見られるデータ構造があると嬉しい
それが ヒープ 
"""

from heapq import heappop, heappush

que = []
q = int(input())

for i in range(q):
    t, h = map(int, input().split())

    if t == 1:
        #ヒープ que に h を追加する
        heappush(que, h)
    else:
        #ヒープでは、最小値が que[0] にあると考えてよい
        #que が空でなくしかも最小値が h 以下なら繰り返す
        while que and que[0] <= h:
            #ヒープ que の最小値を取り出して削除する
            heappop(que)

    print(len(que))