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
  
  string s;
  
  getline(cin,s);
  vector<string>direction={
    "N", "E", "W", "S", "NE", "NW", "SE", "SW"
  };

    
  if (s.substr(0,2) == direction [4]) {
    cout << "SW" << endl;
  }
  
  else if (s.substr(0,2) == direction [5]) {
    cout << "SE" << endl;
  }
  else if (s.substr(0,2)== direction [6]) {
    cout << "NW" << endl;
  }
  else if (s.substr(0,2) == direction [7]) {
    cout << "NE" << endl;
  }
  else if (s.at(0) == direction [0][0]) {
    cout << "S" << endl;
  }
   else if (s.at(0) == direction [1][0]) {
    cout << "W" << endl;
  }
  else if (s.at(0) == direction [2][0]) {
    cout << "E" << endl;
  }
  else if (s.at(0) == direction [3][0]) {
    cout << "N" << endl;
  }
}
