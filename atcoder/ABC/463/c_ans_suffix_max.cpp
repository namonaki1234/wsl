#include <bits/stdc++.h>
using namespace std;

int main() {
    int N;
    cin >> N;

    vector<pair<int, int>> people; // {L, H} にしておくと扱いやすい

    for (int i = 0; i < N; i++) {
        int H, L;
        cin >> H >> L;
        people.push_back({L, H});
    }

    // L の小さい順に並べる
    sort(people.begin(), people.end());

    vector<int> Ls(N);
    vector<int> suffixMax(N);

    for (int i = 0; i < N; i++) {
        Ls[i] = people[i].first;
    }

    // 後ろから累積最大値を作る
    suffixMax[N - 1] = people[N - 1].second;

    for (int i = N - 2; i >= 0; i--) {
        suffixMax[i] = max(suffixMax[i + 1], people[i].second);
    }

    int Q;
    cin >> Q;

    while (Q--) {
        int T;
        cin >> T;

        // L > T となる最初の位置を探す
        int idx = upper_bound(Ls.begin(), Ls.end(), T) - Ls.begin();

        // 問題で必ず誰か残っているなら、そのまま出力してよい
        cout << suffixMax[idx] << '\n';
    }

    return 0;
}