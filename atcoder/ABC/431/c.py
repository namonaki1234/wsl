import io
import sys

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

n, m, k = map(int, input().split())
h = list(map(int, input().split()))
b = list(map(int, input().split()))

h.sort()
b.sort()

h_copy = h.copy()
b_copy = b.copy()
# print(h,b)

cnt = 0
for hi in h_copy:
    for bi in b_copy:
        if hi <= bi:
            cnt += 1
            if cnt >= k:
                print("Yes")
                exit()
            b_copy.remove(bi)
            break

print("No")
