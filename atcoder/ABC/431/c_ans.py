import io
import sys
from itertools import islice

# 下記に標準入力を記載
_InPUT = """\
12 15 12
748 169 586 329 972 529 432 519 408 587 138 249
656 114 632 299 984 755 404 772 155 506 832 854 353 465 387
"""
sys.stdin = io.StringIO(_InPUT)
# ここからコードを記載
# パーツの組み合わせで倒れないのをK体作ることができるならYes,でなければ No を出力
# nが頭パーツの数、mが体パーツの数
# 頭パーツは体パーツの重さ以下
# 体パーツは頭パーツの重さ以上

N, M, K = map(int, input().split())

H = list(map(int, input().split()))
B = list(map(int, input().split()))

# 軽いほうから並べて
H.sort()
B.sort()

# 頭の先頭 K 個と体の末尾 K 個のうち、軽いほうから見てすべてのペアで体の重さが頭の重さ以上ならば Yes 、そうでなければ No
if all(h <= b for h, b in zip(H[:K], B[-K:])):
    print('Yes')
else:
    print('No')

