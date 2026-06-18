#include <bits/stdc++.h>
using namespace std;

const int INF = 1e9;

// 区間最小値を求めるセグメント木
// query(l, r) で区間 [l, r) の最小値を返す
struct SegmentTree {
    int size;              // 葉の数。N以上の2の累乗にする
    vector<int> tree;      // セグメント木本体

    SegmentTree(int n) {
        size = 1;

        // 葉の数を N 以上の2の累乗にする
        while (size < n) {
            size *= 2;
        }

        // 最小値を求めるので、初期値は十分大きい値にする
        tree.assign(2 * size, INF);
    }

    // pos 番目の値を value に更新する
    void update(int pos, int value) {
        // 葉の位置に移動
        pos += size;

        // 葉を更新
        tree[pos] = value;

        // 親に戻りながら、区間の最小値を更新する
        while (pos > 1) {
            pos /= 2;
            tree[pos] = min(tree[2 * pos], tree[2 * pos + 1]);
        }
    }

    // 区間 [l, r) の最小値を求める
    int query(int l, int r) {
        // 葉の位置に移動
        l += size;
        r += size;

        int res = INF;

        // l と r が交差するまで見る
        while (l < r) {
            // l が右の子なら、その区間は答えに使える
            if (l % 2 == 1) {
                res = min(res, tree[l]);
                l++;
            }

            // r が右の子なら、1つ左の区間が答えに使える
            if (r % 2 == 1) {
                r--;
                res = min(res, tree[r]);
            }

            // 親に移動
            l /= 2;
            r /= 2;
        }

        return res;
    }
};

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    int N;
    cin >> N;

    // Y[x] = 「X座標が x の点の Y座標」
    // 0-indexed で扱う
    vector<int> Y(N);

    for (int i = 0; i < N; i++) {
        int x, y;
        cin >> x >> y;

        x--;
        y--;

        // X座標を添字にして、Y座標を保存する
        // これで Y は X の昇順に並んだ配列になる
        Y[x] = y;
    }

    // セグメント木を作る
    SegmentTree seg(N);

    // Y の値をセグメント木に入れる
    for (int i = 0; i < N; i++) {
        seg.update(i, Y[i]);
    }

    int ans = 0;

    for (int x = 0; x < N; x++) {
        int y = Y[x];

        // 自分より左側、つまり区間 [0, x) の Y の最小値を求める
        int min_y_left;

        if (x == 0) {
            // 一番左の点には左側の点が存在しないので、必ず条件を満たす
            min_y_left = INF;
        } else {
            min_y_left = seg.query(0, x);
        }

        // 左側のどの点よりも今の点の Y が小さいなら、
        // 長方形の内部に点が存在しない
        if (y < min_y_left) {
            ans++;
        }
    }

    cout << ans << '\n';

    return 0;
}