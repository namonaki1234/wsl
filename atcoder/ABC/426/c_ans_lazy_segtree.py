import io
import sys

# 下記に標準入力を記載
_InPUT = """\
8 5
2 6
3 5
1 7
5 7
7 8

"""
sys.stdin = io.StringIO(_InPUT)
# ここからコードを記載

INF_NOOP = None  # 遅延の「何もしない」印

class LazySegTree:
    """
    区間代入(assign) と 区間和(sum) を持つ遅延セグ木。
    各葉は「値 val」と「長さ len=1」を持ち、内部節点は和を持つ。
    遅延は「この区間を val に一括代入する」という値（None は NO-OP）。
    """
    def __init__(self, init_vals):  # init_vals: [(val, len=1), ...]
        self.n = len(init_vals)
        self.size = 1
        while self.size < self.n:
            self.size <<= 1
        # 値（区間和）と区間長
        self.val = [0] * (2 * self.size)
        self.len = [0] * (2 * self.size)
        # 遅延（区間一括代入の値／None は未設定）
        self.lazy = [INF_NOOP] * (2 * self.size)

        # 葉をセット
        for i, (v, L) in enumerate(init_vals):
            self.val[self.size + i] = v
            self.len[self.size + i] = L
        for i in range(self.size + self.n, 2 * self.size):
            self.len[i] = 1  # 使わない葉でも長さは 1 にしておく（影響なし）

        # 内部節点を構築（和）
        for i in range(self.size - 1, 0, -1):
            self.val[i] = self.val[i << 1] + self.val[i << 1 | 1]
            self.len[i] = self.len[i << 1] + self.len[i << 1 | 1]

    # 遅延を現在ノードに適用：この区間を f に一括代入
    def _apply(self, k, f):
        if f is INF_NOOP:
            return
        self.val[k] = f * self.len[k]
        self.lazy[k] = f  # composition: 新しい代入が古い代入を上書き

    # 子へ遅延を伝播
    def _push(self, k):
        f = self.lazy[k]
        if f is INF_NOOP:
            return
        self._apply(k << 1, f)
        self._apply(k << 1 | 1, f)
        self.lazy[k] = INF_NOOP

    # 区間 [l, r) を値 f に一括代入
    def range_assign(self, l, r, f):
        def _range_assign(k, nl, nr):
            if r <= nl or nr <= l:
                return
            if l <= nl and nr <= r:
                self._apply(k, f)
                return
            self._push(k)
            mid = (nl + nr) // 2
            _range_assign(k << 1, nl, mid)
            _range_assign(k << 1 | 1, mid, nr)
            self.val[k] = self.val[k << 1] + self.val[k << 1 | 1]

        _range_assign(1, 0, self.size)

    # 区間 [l, r) の和
    def range_sum(self, l, r):
        def _range_sum(k, nl, nr):
            if r <= nl or nr <= l:
                return 0
            if l <= nl and nr <= r:
                return self.val[k]
            self._push(k)
            mid = (nl + nr) // 2
            return _range_sum(k << 1, nl, mid) + _range_sum(k << 1 | 1, mid, nr)

        return _range_sum(1, 0, self.size)

    # 点取得
    def point_get(self, i):
        i += self.size
        # 根から葉へ降りながら遅延を適用
        path = []
        k = i
        while k > 1:
            k >>= 1
            path.append(k)
        for k in reversed(path):
            self._push(k)
        return self.val[i]

    # 点代入（葉を直接セット）
    def point_set(self, i, new_val):
        # 根から葉へ遅延を落としてからセット
        idx = i + self.size
        path = []
        k = idx
        while k > 1:
            k >>= 1
            path.append(k)
        for k in reversed(path):
            self._push(k)

        self.val[idx] = new_val
        # 親を更新
        idx >>= 1
        while idx:
            self.val[idx] = self.val[idx << 1] + self.val[idx << 1 | 1]
            idx >>= 1


n, q = map(int,input().split())

# C++ 版と同様、サイズ n+1 の配列を作る。
# v[0] = {val=0, len=1}、v[1..n] = {val=1, len=1}
init = [(0, 1)] + [(1, 1) for _ in range(n)]
seg = LazySegTree(init)

out = []
for _ in range(q):
    x, y = map(int,input().split())
    # res = sum_{i=0..x} val[i]
    res = seg.range_sum(0, x + 1)
    # [0, x] を 0 に区間代入（対象を全消去）
    seg.range_assign(0, x + 1, 0)
    # 位置 y の現在値に res を加えた値で点代入
    cur_y = seg.point_get(y)
    seg.point_set(y, cur_y + res)
    out.append(str(res))

sys.stdout.write("\n".join(out))
