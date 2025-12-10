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

struct RandomNumberGenerator{
    RandomNumberGenerator() : gen(std::random_device{}()){}
    ll generate(ll l,ll r){
        uniform_int_distribution<ll> dis(l, r);
        return dis(gen);
    }
    mt19937 gen;
} gen;

void solve() {
    int n, m; cin >> n >> m;
    vector<vector<pii>> g(n + 1);
    for (int i = 1;i<m+1;i++) {
        int u, v; cin >> u >> v;
        g[u].push_back(pii(v, i));
        g[v].push_back(pii(u, i));
    }

    vector<bool> bad(m + 1, 0);
    vector<bool> vis(n + 1, 0);
    auto dfs = [&](auto&&self, int u) -> void {
        vis[u] = 1;
        for (auto [v, w] : g[u]) {
            if (vis[v]) continue;
            bad[w] = 1;
            self(self, v);
        }
    };
    dfs(dfs, 1);

    vector<ll> val(m + 1, 0);
    for (int i = 1;i<m+1;i++) {
        if (!bad[i]) {
            val[i] = gen.generate(0, LLONG_MAX);
        }
    }

    
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
    cin >> t;

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
