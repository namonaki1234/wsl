// g++ -D LOCAL c_ans.cpp
// ./a.out

#include <bits/stdc++.h>
using namespace std;

int main(void){
// ==========================================
    // 提出時はここから下のブロックを削除・コメントアウトする
    // ローカル環境でのみ実行されるブロック
#ifdef LOCAL
    string test_input = R"(
5
2 4 5 1 3
4 1 5 2 3
)";
    stringstream ss(test_input);
    cin.rdbuf(ss.rdbuf()); // cin の読み込み先を test_input に変更
#endif
// ==========================================

    int n;
    cin >> n;

    vector<int> a(n);
    for (int i = 0; i < n; i++) {
        cin >> a[i];
    }

    vector<int> b(n);
    for (int i = 0; i < n; i++) {
        cin >> b[i];
    }

    int cnt = 0;
    for (int i = 0; i < n; i++) {
        if ((i + 1) == a[b[i] - 1]) {
            cnt++;
        }
    }

    if (cnt == n) {
        cout << "Yes" << endl;
    } else {
        cout << "No" << endl;
    }

    return 0;
}