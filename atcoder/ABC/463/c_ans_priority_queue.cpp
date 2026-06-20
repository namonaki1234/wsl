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

    // L が大きい順
    sort(people.rbegin(), people.rend());

    int Q;
    cin >> Q;

    vector<pair<int, int>> queries; // {T, query_index}

    for (int i = 0; i < Q; i++) {
        int T;
        cin >> T;
        queries.push_back({T, i});
    }

    // T が大きい順
    sort(queries.rbegin(), queries.rend());

    priority_queue<int> pq;
    vector<int> ans(Q);

    int p = 0;

    for (auto [T, qi] : queries) {
        // L > T を満たす人を追加
        while (p < N && people[p].first > T) {
            int H = people[p].second;
            pq.push(H);
            p++;
        }

        ans[qi] = pq.top();
    }

    for (int i = 0; i < Q; i++) {
        cout << ans[i] << '\n';
    }

    return 0;
}