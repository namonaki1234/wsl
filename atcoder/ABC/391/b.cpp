#include <bits/stdc++.h>
using namespace std;
// #include <atcoder/all>
// using namespace atcoder;
#define rep(i, n) for (int i = 0; i < (int)(n); i++)
#define rep1(i, n) for (int i = 1; i < (int)(n+1); i++)
using ll = long long;
using P = pair<int,int>;
using Graph = vector<vector<int>>;
// using mint = modint1000000007;

int main() {
    int S, T;
    
    cin >> S >> T;
    cin.ignore(); // 改行を消去

    vector<string> grid_S(S);
    for (int i = 0; i < S; ++i) {
        getline(cin, grid_S[i]);
    }
     vector<string> grid_T(T);
    for (int i = 0; i < T; ++i) {
        getline(cin, grid_T[i]);
    }
    int count = 0;

    for (int i = 0; i < S - 1; ++i) {
        for (int j = 0; j < S - 1; ++j) {
            if (grid_S[i][j] != grid_T[i][j] ){
              continue;
            }
            else if (grid_S[i][j] == grid_T[i][j] ){
              count++;

              if (count == T*T){
                int a = i;
                int b = j;
                cout <<  a << " " << b << endl;
                break;
              }
             
            
                
            }
        }
    }

    return 0;
}
