import io
import sys
from collections import deque
import itertools
# from atcoder.dsu import DSU
# 再帰の上限を増やす（Pythonでは必須のおまじない）
sys.setrecursionlimit(10**6)

# 下記に標準入力を記載
_InPUT = """\
10 2
3 3 4 1 1 3 3 1 5 1
"""
sys.stdin = io.StringIO(_InPUT)
# ここからコードを記載
#n種類のアイテム、はじめ高橋君はアイテム1を持っている
#高橋にはM人と友達がいる。i人目の友達にアイテムA_iを渡すと、アイテムB_iをもらうことができる
#高橋が手に入れることのできるアイテムはアイテム1を含めて何種類

n,k = map(int, input().split())
a = [int(x) for x in input().split(" ")]

for i in range(k):
    max_a = max(a)
    for j in range(a.count(max_a)):
        a[a.index(max_a)] = 0
    # a.remove(max_a)
    print(a)
    
print(sum(a))