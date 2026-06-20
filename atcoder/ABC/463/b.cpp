#include <iostream>
#include <vector>

using namespace std;
// Aのasciiは65、Cは67だからA-Eは65-69
int main() {
    int n;
    char x;
    cin >> n >> x;

    for (int i = 0; i < n; ++i) {
        string s;
        cin >> s;
        int index = int(x) - 65;
        if (s[index] == 'o'){
            cout << "Yes" << endl;
            return 0;
        }
    }
    cout << "No" << endl;

    return 0;
}