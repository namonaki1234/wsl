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
5 3 3
1 30
1 40
1 50
2 10
3 20
)";
    stringstream ss(test_input);
    cin.rdbuf(ss.rdbuf()); // cin の読み込み先を test_input に変更
#endif
// ==========================================

    int n, k, m;
    cin >> n >> k >> m;

    // t[color] に、その色の宝石の価値をまとめる
    // Pythonの t = [[] for _ in range(N + 1)] に相当
    vector<vector<long long>> t(n + 1);

    for (int i = 0; i < n; i++) {
        int c;
        long long v;
        cin >> c >> v;
        t[c].push_back(v);
    }

    // top  : 各色の代表宝石
    // tail : 代表以外の宝石
    vector<long long> top;
    vector<long long> tail;

    // 色ごとに処理
    for (int i = 1; i <= n; i++) {
        // その色の宝石が1つもないならスキップ
        if (t[i].empty()) {
            continue;
        }

        // 価値が高い順（降順）に並べる
        sort(t[i].rbegin(), t[i].rend());

        // 一番価値が高いものをその色の代表にする
        top.push_back(t[i][0]);

        // 残りは代表以外の候補に入れる (Pythonの tail += values[1:])
        for (size_t j = 1; j < t[i].size(); j++) {
            tail.push_back(t[i][j]);
        }
    }

    // 代表宝石を価値が高い順に並べる
    sort(top.rbegin(), top.rend());

    long long answer = 0;

    // 代表のうち、上位M個は必ず選ぶ
    for (int i = 0; i < m; i++) {
        answer += top[i];
    }

    // M個に入らなかった代表は、普通の候補(tail)に戻す
    for (size_t i = m; i < top.size(); i++) {
        tail.push_back(top[i]);
    }

    // 残り候補を価値が高い順に並べる
    sort(tail.rbegin(), tail.rend());

    // 代表M個を選んだので、あとK-M個をtailから選ぶ
    for (int i = 0; i < k - m; i++) {
        answer += tail[i];
    }

    cout << answer << endl;

    return 0;
}