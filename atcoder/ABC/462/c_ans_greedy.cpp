#include <bits/stdc++.h>
using namespace std;

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    int N;
    cin >> N;

    vector<int> Y(N);

    for (int i = 0; i < N; i++) {
        int x, y;
        cin >> x >> y;
        Y[x - 1] = y;
    }

    int ans = 0;
    int min_y = N + 1;

    for (int y : Y) {
        if (y < min_y) {
            ans++;
        }
        min_y = min(min_y, y);
    }

    cout << ans << '\n';

    return 0;
}