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


import typing


class DSU:
    '''
    Implement (union by size) + (path halving)

    Reference:
    Zvi Galil and Giuseppe F. Italiano,
    Data structures and algorithms for disjoint set union problems
    '''

    def __init__(self, n: int = 0) -> None:
        self._n = n
        self.parent_or_size = [-1] * n

    def merge(self, a: int, b: int) -> int:
        assert 0 <= a < self._n
        assert 0 <= b < self._n

        x = self.leader(a)
        y = self.leader(b)

        if x == y:
            return x

        if -self.parent_or_size[x] < -self.parent_or_size[y]:
            x, y = y, x

        self.parent_or_size[x] += self.parent_or_size[y]
        self.parent_or_size[y] = x

        return x

    def same(self, a: int, b: int) -> bool:
        assert 0 <= a < self._n
        assert 0 <= b < self._n

        return self.leader(a) == self.leader(b)

    def leader(self, a: int) -> int:
        assert 0 <= a < self._n

        parent = self.parent_or_size[a]
        while parent >= 0:
            if self.parent_or_size[parent] < 0:
                return parent
            self.parent_or_size[a], a, parent = (
                self.parent_or_size[parent],
                self.parent_or_size[parent],
                self.parent_or_size[self.parent_or_size[parent]]
            )

        return a

    def size(self, a: int) -> int:
        assert 0 <= a < self._n

        return -self.parent_or_size[self.leader(a)]

    def groups(self) -> typing.List[typing.List[int]]:
        leader_buf = [self.leader(i) for i in range(self._n)]

        result: typing.List[typing.List[int]] = [[] for _ in range(self._n)]
        for i in range(self._n):
            result[leader_buf[i]].append(i)

        return list(filter(lambda r: r, result))
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






for _ in range(n - 1):
    p, c = map(int, input().split())
     # 双方向のペアを追加
    pairs.append((p, c))



q = int(input())
query = []
for i in range(q):
 s = tuple(map(int, input().split()))
 t = tuple(map(int, input().split()))

 query.append((s, t))

# print(query)
