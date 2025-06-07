#include <iostream>
#include <vector>
#include <algorithm>

int main() {
    using namespace std;
    unsigned H, W;
    cin >> H >> W;
    vector A(H, vector<unsigned long>(W));
    for (auto&& row : A)
        for (auto&& a : row)
            cin >> a;

    // HW ビットの整数でドミノが置かれたマスを管理する
    vector<unsigned> possible_domino{0}; // はじめは何も置かれていない盤面だけ

    const unsigned domino_vertical = (1U << W) + 1, domino_horizontal = 3;

    for (unsigned i = 0, bit = 0; i < H; ++i)
        for (unsigned j = 0; j < W; ++j, ++bit) { // 各マスについて
            vector<unsigned> tmp;
            for (const auto b : possible_domino) {
                // 縦に置けるなら置いたものを追加する
                if (j + 1 < W && !(b & domino_horizontal << bit))
                    tmp.emplace_back(b | domino_horizontal << bit);
                // 横に置けるなら置いたものを追加する
                if (i + 1 < H && !(b & domino_vertical << bit))
                    tmp.emplace_back(b | domino_vertical << bit);
            }
            ranges::move(tmp, back_inserter(possible_domino));
        }

    unsigned long ans = 0;
    for (const auto b : possible_domino) { // 各ドミノの置き方について
        unsigned long now = 0;
        for (unsigned i = 0, bit = 0; i < H; ++i)
            for (unsigned j = 0; j < W; ++j, ++bit)
                if (1 & ~b >> bit) // ドミノが置かれてないマスに
                    now ^= A[i][j]; // 書かれた整数の xor を求める
        ans = max(ans, now); // 大きければ更新
    }

    cout << ans << endl;
    return 0;
}
