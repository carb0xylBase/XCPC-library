#include <bits/stdc++.h>
using namespace std;
#define DEBUG 1
using ll = long long;
using pii = pair<ll,ll>;
const ll N = 2000000;
const ll INF = 5e18;
const ll MOD = 1e9 + 7;

struct Kmp {
    // a 为文本串, b 为模式串
    vector<int> get_pos(string& a,string& b) {
        vector<int> res;
        string z = b;
        z.push_back('@'); // 隔断符
        z += a;

        vector<int> pi(z.size(),0);
        for (int i = 1;i<z.size();i++) {
            int j = pi[i - 1];
            while (j > 0 && z[j] != z[i]) {
                j = pi[j-1];
            }
            if (z[j] == z[i]) j++;
            pi[i] = j;
        }

        for (int i = b.length() + 1;i<z.length();i++) {
            if (pi[i] >= b.length()) {
                res.push_back(i - b.length() * 2);
            }
        }

        return res;
    }
};