import io
import sys
from collections import defaultdict,deque,Counter
import numpy as np
from itertools import combinations

# 下記に標準入力を記載
_InPUT = """\
999999999
"""
sys.stdin = io.StringIO(_InPUT)
# ここからコードを記載
sys.setrecursionlimit(1000000)

# https://atcoder.jp/contests/abc114/submissions/20050058

N = int(input())

ans = 0

# 数 n について調べ、末尾に数字を追加した数を再帰的に調べる関数。
# use3, use5, use7 はそれぞれ 3, 5, 7 を使ったかというフラグ。
def dfs(n, use3, use5, use7):
    global ans
    # N を超えていたら打ち切って終了する。
    if n > N:
        return
    # 3 種類全てを使っていたら答えに加算する。
    if use3 and use5 and use7:
        ans += 1
    # 末尾に 3, 5, 7 を付けた数を再帰的に調べる。
    dfs(10*n+3, True, use5, use7)
    dfs(10*n+5, use3, True, use7)
    dfs(10*n+7, use3, use5, True)

# 何もない状態（値としては 0）から呼び出す。
dfs(0, False, False, False)

# 答えを出力する。
print(ans)

'''
いつ終了するか？
それぞれの dfs(n, ...) の中で最初にこうチェックしますよね：

python
コピーする
編集する
if n > N:
    return
たとえば n = 3333333 などになったときに、N = 1000 を超えていたらここで終了：

python
コピーする
編集する
if n > 1000:  # True
    return    # ← ここで関数から抜ける（＝終了）
🔁 この「return」で戻るということは…
dfs(3333333) が return → 呼び出し元 dfs(333333) に戻る

dfs(333333) の処理の残りに移る（次の dfs(10*n+5) の行など）

それがまた N > 1000 で return → 上の dfs(33333) に戻る

というふうに「スタックが1つずつ解かれて戻っていく」のです。
'''