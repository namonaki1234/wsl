#include <iostream>
#include <vector>

using namespace std;

int main() {
    // コンテスト向けの入出力高速化
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    int n;
    if (!(cin >> n)) return 0;

    // Pythonの ans = [[] for _ in range(n+1)] に相当
    // ans[i] には、人iにギフトを送った人の番号が格納される
    vector<vector<int>> ans(n + 1);

    for (int i = 1; i <= n; ++i) {
        int k;
        cin >> k; // 最初に送った人数 k を受け取る
        
        for (int j = 0; j < k; ++j) {
            int a_j;
            cin >> a_j; // 受け取った人の番号を1つずつ受け取る
            
            // 受け取った人(a_j)のリストに、送った人(i)を追加
            ans[a_j].push_back(i);
        }
    }

    // 1番の人から n番の人まで順に出力
    // Pythonの for j in range(1, len(ans)): と同じ範囲
    for (int i = 1; i <= n; ++i) {
        // まず送ってくれた人数 (リストの長さ) を出力
        cout << ans[i].size();
        
        // 続いて、送ってくれた人の番号を空白区切りで出力
        // Pythonの *ans[i] に相当する処理
        for (int sender : ans[i]) {
            cout << " " << sender;
        }
        cout << "\n";
    }

    return 0;
}