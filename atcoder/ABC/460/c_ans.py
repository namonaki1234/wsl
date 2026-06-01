import io
import sys
from collections import Counter, deque
import itertools
import bisect
import re
# from atcoder.dsu import DSU

# 再帰の上限を増やす（Pythonでは必須のおまじない）
sys.setrecursionlimit(10**6)

# 下記に標準入力を記載
_InPUT = """\
4 5
4 2 1 8
14 9 3 2 9
"""
sys.stdin = io.StringIO(_InPUT)

# ここからコードを記載
# n個のシャリ、m個のネタ
# i番目のシャリの重さa_i,j番目のネタの重さb_j
# ネタの重さ<=シャリの重さ*2

n,m = map(int,input().split())

a = [int(x) for x in input().split()]
b = [int(x) for x in input().split()]

# 【ポイント1】まずは両方とも小さい順（昇順）にソートする
a.sort()
b.sort()

# 現在見ているシャリのインデックス(i) と ネタのインデックス(j)
i = 0
j = 0
ans = 0

# 【ポイント2】一筆書きでリストを左から右へ1回だけ見る
while i < n and j < m:
    # 条件：ネタの重さがシャリの2倍以下なら、お寿司を1つ作る
    if b[j] <= 2 * a[i]:
        ans += 1  # お寿司が1つ完成
        i += 1    # 次のシャリへ
        j += 1    # 次のネタへ
    else:
        # ネタ(b[j])が重すぎる場合、このシャリ(a[i])ではどう頑張っても乗せられない。
        # ネタはキープしたまま、シャリだけを「次の少し大きいシャリ」に進めてみる。
        i += 1

print(ans)