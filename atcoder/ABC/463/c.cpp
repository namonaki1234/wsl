#include <bits/stdc++.h>
#include <iostream>
#include <vector>
#include <algorithm> // sort() を使うために必要

using namespace std;

int main() {
    int n;
    cin >> n;

    vector<pair<int,int>> hl;
    // map<int, int, greater<int>> hl_map{};
    // vector<int> hl{};

    for (int i = 0; i < n; i++) {
        int h, l;
        cin >> h >> l;
        hl.push_back({h, l});
    }

    int q;
    cin >> q;
    vector<int> t(q);
    for (int i = 0; i < q; ++i) cin >> t[i];

    // for (auto [key, val]: hl_map) cout << "[" << key << "]=" << val << endl;
    // for (auto [h,l]: hl) cout << h << l << " ";
    // for (auto it: t) cout << it << " ";
    sort(hl.rbegin(), hl.rend());
    int ans = 0;
    for (int i = 0;i < q; ++i) {
        // cout << hl[i].first << t[i]+0.5 << endl;
        if  (hl[i].first >= (t[i] + 0.5)) {
            // ans = max(ans, hl[i].first);
            cout << hl[i].first << endl;
        }
        cout << ans << endl;
    }
    return 0;
}