#include <bits/stdc++.h>
using namespace std;

int main() {
    int N;
    cin >> N;

    vector<pair<int, int>> cand;

    for (int i = 0; i < N; i++) {
        int H, L;
        cin >> H >> L;

        // 今の人より背が低く、かつ今の人より早く退室する候補は不要
        while (!cand.empty() && cand.back().first <= H) {
            cand.pop_back();
        }

        cand.push_back({H, L});
    }

    int Q;
    cin >> Q;

    for (int i = 0; i < Q; i++) {
        int T;
        cin >> T;

        // cand は L 昇順になっている
        // L > T となる最初の候補を探す
        int left = 0;
        int right = cand.size();

        while (left < right) {
            int mid = (left + right) / 2;

            if (cand[mid].second > T) {
                right = mid;
            } else {
                left = mid + 1;
            }
        }

        cout << cand[left].first << endl;
    }

    return 0;
}