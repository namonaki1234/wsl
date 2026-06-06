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
4 5
)";
    stringstream ss(test_input);
    cin.rdbuf(ss.rdbuf()); // cin の読み込み先を test_input に変更
#endif
// ==========================================

    int a, d;
    cin >> a >> d;

    if (a <= d) {
        cout << "Yes" << endl;
    } else {
        cout << "No" << endl;
    }

    return 0;
}