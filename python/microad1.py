# 基本問題仕様

# プログラムの標準入力に、N-1個のidのペアが与えられます。
# その後クエリとしてQ個のidのペアが与えられるので、
# 各ペアが最初に与えられたN-1個のペアに含まれていた場合は `Yes` を、
# 含まれていなかった場合は `No` を標準出力に出力してください。
# **この時、クエリで与えられたidの順序を区別せずに判定を行ってください。**

## 入力形式

# 入力は、下記の形式でプログラムの標準入力から与えられます。

# ```
# N
# p_1, c_1
# p_2, c_2
# ...
# p_{N-1}, c_{N-1}
# Q
# u_1, v_1
# u_2, v_2
# ...
# u_Q, v_Q
# ```

# * N: ペアに出現するidの種類数
# * p_i, c_i: idのペア
#   * `1 <= p_i <= N`
#   * `1 <= c_i <= N`
# * Q: 与えられるクエリの数
#   * `1 <= Q <= N*(N-1)`
# * u_i, v_i: 判定対象のidのペア
#   * `1 <= u_i <= N`
#   * `1 <= v_i <= N`

# ## 出力

# 各クエリに対する判定結果、`Yes` or `No` を標準出力に表示してください。

import io
import sys
from collections import defaultdict,deque,Counter
import numpy as np
from itertools import combinations

# 下記に標準入力を記載
_InPUT = """\
3
1 2
2 3
3
1 2
1 3
2 1

"""
sys.stdin = io.StringIO(_InPUT)
# ここからコードを記載

n = int(input())
pairs = set()
for _ in range(n - 1):
    p, c = map(int, input().split())
     # 双方向のペアを追加
    pairs.add((p, c))
    pairs.add((c, p))

query = [tuple(map(int, input().split())) for _ in range(int(input()))]
# print(pairs)
print("\n".join("Yes" if (u, v) in pairs else "No" for u, v in query))