// g++ -D LOCAL c_ans.cpp
// ./a.out

#include <iostream>
#include <string>

int main() {
    std::string s;
    // 標準入力から文字列を受け取る
    std::cin >> s;

    std::string ans = "";
    
    // 範囲for文で文字列sから1文字ずつ取り出して処理
    for (char c : s) {
        // 文字 c が 'a' から 'z' の間に含まれていない場合
        if (c < 'a' || c > 'z') {
            ans += c;
        }
    }

    // 結果を出力
    std::cout << ans << "\n";

    return 0;
}