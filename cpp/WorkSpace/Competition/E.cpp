#include <bits/stdc++.h>
using namespace std;
typedef long long ll;
#define DEBUG 1
const ll N = 2000000;
const ll MOD = 998244353;
const ll INF = 1e18; 

string s;
bool check(string t) {
    ll L = 0,R = 0;
    char opt = '+';
    int j = 0;
    while (1) {
        int i;
        ll res = 0;
        for (i = j;;i++) {
            if (t[i] >= '0' && t[i] <= '9') {
                res *= 10;res += t[i] - '0';
            } else {
                if (opt == '+') {

                } else {
                    res *= -1;
                }
                opt = t[i];
                break;
            }
        }

        if (i == j) {
            return 0;
        }

        if (i - j > 10) {
            return 0;
        }

        if (i - j > 1 && t[j] == '0') {
            return 0;
        }

        L += res;

        if (opt == '=') {
            break;
        }
        j = i + 1;
    }

    int cur = t.find('=') + 1;
    if (cur >= t.length()) return 0;
    opt = '+';
    while (1) {
        int i;
        ll res = 0;
        for (i = cur;i<t.length();i++) {
            if (t[i] >= '0' && t[i] <= '9') {
                res *= 10;res += t[i] - '0';
                if (i == t.length() - 1) {
                    if (opt == '-') {
                        res *= -1;
                    }
                }
            } else {
                if (opt == '+') {

                } else {
                    res *= -1;
                }

                opt = t[i];
                break;
            }
        }

        if (cur == i) {
            return 0;
        }

        if (i - cur > 1 && t[cur] =='0') {
            return 0;
        }

        if (i - cur > 10) {
            return 0;
        }

        R += res;

        cur = i + 1;
        if (cur >= t.length()) {
            break;
        }
    }

    if (L == R) {
        return 1;
    }
    return 0;
}

void solve() {
    cin >> s;
    if (check(s)) {
        cout << "Correct" << endl;
        return;
    }

    for (int i = 0;i<s.length();i++) {
        if (s[i] >= '0' && s[i] <= '9') {
            string t = s;
            t.erase(i,1);
            for (int j = 0;j<t.length();j++) {
                string z = t;
                z.insert(z.begin()+j,s[i]);
                // cout << z << endl;
                if (check(z)) {
                    cout << z << endl;
                    return;
                }
            }
        }
    }

    cout << "Impossible" << endl;
    return;
}

signed main() {
    // freopen("input.txt","r",stdin);
    ios::sync_with_stdio(0),cin.tie(0),cout.tie(0);
    int _ = 1;
    // cin >> _;
    while (_--){
        solve();
    }
    return 0;
}
