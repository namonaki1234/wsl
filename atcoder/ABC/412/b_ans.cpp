#include <bits/stdc++.h>
using namespace std;

int main() {
	string s, t;
	cin >> s >> t;
	bool res = true;
	int n = ssize(s);
	for (int i = 1; i < n; i++) {
		if (isupper(s[i])) {
			if (t.find(s[i - 1]) == string::npos) res = false;
		}
	}
	cout << (res ? "Yes" : "No") << '\n';
}
