#define DEBUG 1
#define FUCK cout << "fuck" << endl;
#if DEBUG
    #include "all.hpp"
#else
    #include <bits/stdc++.h>
#endif

using namespace std;
using ll = long long;
using pii = pair<int,int>;
using pll = pair<ll,ll>;
using db = long double;
using pdd = pair<db, db>;

const ll N = 50000000;
const ll INF = 5e18;
const ll MOD = 1e9 + 7;

namespace Z {
    string s; 
    int z[N];
    void init() {
        z[0] = 0;
        int n = s.size();
        for (int i = 1,l = 0,r = 0;i<n;i++) {
            if (i <= r && z[i - l] < r - i + 1) {
                z[i] = z[i - l];
            } else {
                z[i] = max(0, r - i + 1);
                while (i + z[i] < n && s[z[i]] == s[i + z[i]]) ++ z[i];
            }
            if (i + z[i] - 1 > r) l = i, r = i + z[i] - 1;
        }
        return;
    }
};
using namespace Z;

void solve() {
    string a, b;
    cin >> a >> b;
    s = b;
    init();
    ll ans = 0;
    z[0] = b.size();
    for (int i = 0;i<b.size();i++) {
        ans ^= 1ll * (z[i] + 1) * (i + 1);
    }
    cout << ans << endl;
    ans = 0;
    s = b + "#" + a;
    init();
    vector<int> p(a.size(), 0);
    for (int i = b.size()+1;i<s.size();i++) {
        ans ^= 1ll * (z[i] + 1) * (i - b.size());
    }
    cout << ans << endl;
    return;
}

signed main() {
#if DEBUG
    freopen("input.txt", "r", stdin);
    auto start_time = chrono::steady_clock::now();
#else
    ios::sync_with_stdio(false);
#endif
    cin.tie(nullptr);

    int t = 1;
    // cin >> t;

    while (t--) {
        solve();
    }

#if DEBUG
    auto end_time = chrono::steady_clock::now();
    auto diff = chrono::duration_cast<chrono::milliseconds>(end_time - start_time);
    cerr << "Time: " << diff.count() << " ms" << endl;

#endif

    return 0;
}
