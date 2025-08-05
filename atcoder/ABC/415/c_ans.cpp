#include<bits/stdc++.h>

using namespace std;

int main(){
  int t;
  cin >> t;
  while(t>0){
    t--;
    int n;
    cin >> n;
    string s;
    cin >> s;
    s="0"+s;
    vector<int> ok(1<<n,0);
    ok[0]=1;
    for(int i=0;i<(1<<n);i++){
      if(ok[i]==0){continue;}
      for(int j=0;j<n;j++){
        if(i&(1<<j)){continue;}
        int next=(i|(1<<j));
        if(s[next]=='0'){ ok[next]=1; }
      }
    }
    if(ok[(1<<n)-1]){cout << "Yes\n";}
    else{cout << "No\n";}
  }
  return 0;
}
