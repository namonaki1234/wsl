# # 応用問題の背景

# 広告を配信する際には必ず広告を表示する先=広告枠があります。
# 広告枠はニュースサイトやスマートフォンアプリなど様々なメディアに存在し多様な種類や属性が存在します。

# ![広告配信面のイメージ](img/ad_delivery.png)

# メディアの情報を収集して属性（トピック）について調査することをトピック分析と言います。
# トピックによって広告の配信効果が変わる場合があるため、マイクロアドでも配信を切り替えするためのメディアのトピック情報を管理しています。

# マイクロアドで管理しているトピックは階層構造になっています。
# たとえばハイブリッドの記事だと
# `クルマ＆乗り物/クルマ/ハイブリッド`

# トラックなどの記事だと
# `クルマ＆乗り物/クルマ/トラック`

# 航空会社のマイル関連のサイトだと下記のようなトピック階層になっているかもしれません。
# `クルマ＆乗り物/飛行機/マイル`

# 広告を配信する際にはハイブリッドを取り扱ったWebページに配信したい場合もありますが、クルマというハイブリッドに比べて広いトピックに対して配信したい場合があるかもしれません。
# それを実現するために広告配信システムでは配信対象のトピックの親子関係を解決して配信するべきかどうか高速に判断する必要があります。

# 今回コーディング問題としてこの課題を取り上げました。

# # 応用問題１仕様

# プログラムの標準入力に、
# - N個のカテゴリ名
# - N-1個のカテゴリidの親子関係

# からなる、全てのカテゴリがグラフ上で連結となるカテゴリ木が与えられます。

# その後クエリとしてQ個のカテゴリidが与えられるので、
# 「**根ノードを除く**親カテゴリ名」と「与えられたidのカテゴリ名」を `/` 区切りで結合し、標準出力に出力してください。
# (より具体的な出力形式は、問題入出力例を参照してください)

# ## 入力形式

# 入力は、下記の形式でプログラムの標準入力から与えられます。

# ```
# N
# A_1, A_2, ..., A_N
# p_1, c_1
# p_2, c_2
# ...
# p_{N-1}, c_{N-1}
# Q
# u_1
# u_2
# ...
# u_Q
# ```

# * N: 事前に与えられるカテゴリの数
# * A_i: `カテゴリid=i` のカテゴリ名(A_1は根ノード)
# * p_i, c_i: 親子関係を持つ根ノード・カテゴリidのペア(pが親、cが子の関係)
#   * `1 <= p_i <= N`
#   * `2 <= c_i <= N`
# * Q: 与えられるクエリの数
# * u_i: カテゴリid
#   * `2 <= u_i <= N`

# ## 出力

# 各クエリで与えられたカテゴリidのカテゴリ名(根ノード直前の親カテゴリ含む)を標準出力に表示してください。

import io
import sys
from collections import defaultdict,deque,Counter
import numpy as np
from itertools import combinations

# 下記に標準入力を記載
_InPUT = """\
5
root エンタメ 映画 アニメ 邦画
1 2
2 3
2 4
3 5
2
2 3
3 4
4 5
2 3
"""
sys.stdin = io.StringIO(_InPUT)
# ここからコードを記載

n = int(input())
a = list(input().split())
# print(a)
# num = range(1, n + 1)
pairs = list()

N,M = len(a), n - 1

root = [i for i in range(N+1)] # 要素の根を格納する配列
rank = [0]*(N+1) # 木のランク（高さ）を格納する配列
size = [1]*(N+1) # 木のサイズ（要素数）を格納する配列

def union(x,y):
  root_x = find(x)
  root_y = find(y)
  if root_x != root_y:
    if rank[root_x] >= rank[root_y]:
      if rank[root_x] == rank[root_y]:
        rank[root_x] += 1
      root[root_y] = root_x
      size[root_x] += size[root_y]
    else:
      root[root_x] = root_y
      size[root_y] += size[root_x]

def find(x):
  tmp_root = x
  while True:
    if tmp_root == root[tmp_root]:
      break
    tmp_root = root[tmp_root]
  root[x] = tmp_root
  return tmp_root

# find(A) == find(B) ならばAとBは同じグループに属する
# size[find(A)] でAを含むグループの要素数を得られる




for _ in range(n - 1):
    p, c = map(int, input().split())
     # 双方向のペアを追加
    pairs.append((p, c))

for i in range(M):
  union(pairs[i][0], pairs[i][1])


q = int(input())
query = []
for i in range(q):
 s = tuple(map(int, input().split()))
 t = tuple(map(int, input().split()))

 query.append((s, t))

# print(query)
