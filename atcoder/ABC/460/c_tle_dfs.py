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



# ネタがすでに使われたかどうかを管理するリスト（すべてFalseで初期化）
used_b = [False] * m

def dfs(shari_idx):
    """
    shari_idx 番目以降のシャリを使って作れるお寿司の最大数を返す関数
    """
    # ベースケース：すべてのシャリを調べ終わったら終了（これ以上作れないので0を返す）
    if shari_idx == n:
        return 0
    
    # パターン1：現在のシャリ (a[shari_idx]) を【使わない（スキップする）】場合
    max_sushi = dfs(shari_idx + 1)  
    
    # パターン2：現在のシャリと、【条件を満たすまだ使っていないネタ】を組み合わせる場合
    for j in range(m):
        if not used_b[j] and b[j] <= 2 * a[shari_idx]:
            
            # --- 探索に進む前の準備 ---
            # ネタjを「使用済み」にする
            used_b[j] = True
            
            # --- 深さ優先探索 ---
            # お寿司が1つ完成したので 1 を足し、次のシャリ (shari_idx + 1) の探索へ進む
            current_sushi = 1 + dfs(shari_idx + 1)
            
            # 最大値を更新
            max_sushi = max(max_sushi, current_sushi)
            
            # --- バックトラック（状態を戻す） ---
            # 別の組み合わせ（他のネタを選ぶパターン）も試すために、
            # 探索から戻ってきたらネタjを「未使用」に戻す
            used_b[j] = False
            
    return max_sushi

# 最初のシャリ(index 0)から探索をスタートし、結果を出力
ans = dfs(0)
print(ans)