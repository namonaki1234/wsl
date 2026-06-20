#include <bits/stdc++.h>
using namespace std;

int main() {
    int N;
    cin >> N;

    vector<pair<int, int>> people; // {L, H}

    for (int i = 0; i < N; i++) {
        int H, L;
        cin >> H >> L;
        people.push_back({L, H});
    }

    sort(people.begin(), people.end()); // L 昇順

    int Q;
    cin >> Q;

    vector<pair<int, int>> queries; // {T, query_index}

    for (int i = 0; i < Q; i++) {
        int T;
        cin >> T;
        queries.push_back({T, i});
    }

    sort(queries.begin(), queries.end()); // T 昇順

    multiset<int> heights;

    // 最初は全員いる
    for (auto [L, H] : people) {
        heights.insert(H);
    }

    vector<int> ans(Q);

    int p = 0;

    for (auto [T, qi] : queries) {
        // T 分までに退室した人を消す
        while (p < N && people[p].first <= T) {
            int H = people[p].second;

            auto it = heights.find(H);
            heights.erase(it); // 同じ身長が複数いるので、1個だけ消す

            p++;
        }

        ans[qi] = *heights.rbegin();
    }

    for (int i = 0; i < Q; i++) {
        cout << ans[i] << '\n';
    }

    return 0;
}