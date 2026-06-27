#include <iostream>
#include <vector>
#include <algorithm>

using namespace std;

int main() {
    int h,w,count=0;
    cin >> h >> w;
    vector<string> c(h);

    for (int i=0;i<h;i++){
        cin >> c[i];
    }

    // 上下左右の#の位置indexを記録する変数を用意
    int top = h-1,bottom=0,left=w-1,right=0;
    // vector<string> c_copy = c;

    for (int i=0;i<h;i++){
        for (int j=0;j<w;j++){
            if (c[i][j] == '#') {
                top = min(top,i);
                bottom = max(bottom,i);
                left = min(left,j);
                right = max(right,j);
            }
        }
    }

    
    for (int i=top;i<=bottom;i++){
        for (int j=left;j<=right;j++){
            cout << c[i][j];
        }
        cout << endl;
    }

    // for (int i=0;i<c.size();i++){
    //     cout << c[i] << endl;
    // }
    
    
    return 0;
}