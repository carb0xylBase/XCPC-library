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
using i128 = __int128_t;

const ll N = 2000000;
const ll INF = 5e18;
const ll MOD = 1e9 + 7;

void solve() {
    int n; cin >> n;
    struct C {
        ll x, y, r;
        int idx;
    };
    vector<C> a(n);
    vector<C> b;
    int cnt_idx = 0;
    for (auto&v : a) {
        cin >> v.x >> v.y >> v.r;
        v.idx = cnt_idx ++;
    }
    b = a;
    vector<int> fa(n, -1);

    set<pll> in, de;
    sort(a.begin(), a.end(), [&](C A, C B) {
        if (A.x - A.r != B.x - B.r) return A.x - A.r < B.x - B.r;
        if (A.y != B.y) return A.y > B.y;
        return A.r > B.r;
    });
    for (auto v : a) {
        cout << v.x << ' ' << v.y << ' ' << v.r << endl;
        while (!de.empty() && de.begin()->first <= v.x - v.r) {
            in.erase(pll(b[de.begin()->second].y + b[de.begin()->second].r, de.begin()->second));
            in.erase(pll(b[de.begin()->second].y, de.begin()->second));
            de.erase(de.begin());
        }
        if (!in.empty()) {
            auto it = in.lower_bound(pll(v.y + v.r, -INF));
            if (it != in.end()) {
                int j = it->second;
                if (b[j].y == it->first) {
                    fa[v.idx] = fa[j];
                } else {
                    fa[v.idx] = j;
                }
            }
        }
        de.insert(pll(v.x + v.r, v.idx));
        in.insert(pll(v.y + v.r, v.idx));
        in.insert(pll(v.y, v.idx));
    }

    vector<vector<int>> g(n);
    for (int i = 0;i<n;i++) {
        if (fa[i] != -1) {
            g[fa[i]].push_back(i);
        }
    }
    ll ans = 0;
    for (int i = 0;i<n;i++) {
        if (fa[i] != -1) continue;
        auto dfs = [&](auto&&self, int u) -> ll {
            ll tot = b[u].r * b[u].r;
            for (auto v : g[u]) {
                tot -= self(self, v);
            }
            return tot;
        };
        ans += dfs(dfs, i);
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
