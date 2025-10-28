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

const ll N = 2000000;
const ll INF = 5e18;
const ll MOD = 998244353;

namespace PAM {
    string s;
    int n, fail[N], len[N], sum[N], ch[N][26], tot;
    int cnt[N], in[N], cnt2[N];

    int getFail(int x, int i) {
        while (i - len[x] - 1 < 0 || s[i - len[x] - 1] != s[i]) {
            x = fail[x];
        }
        return x;
    }

    // 先给 s 赋值
    void init() {
        tot = 1;
        n = s.size();
        fail[0] = 1, len[1] = -1;
        int cur = 0;
        for (int i = 0;i<n;i++) {
            int pos = getFail(cur, i);
            if (!ch[pos][s[i] - '0']) {
                fail[++tot] = ch[getFail(fail[pos], i)][s[i]-'0'];
                ch[pos][s[i] - '0'] = tot;
                len[tot] = len[pos] + 2;
                sum[tot] = sum[fail[tot]] + 1;
            }
            cur = ch[pos][s[i] - '0'];
            cnt[cur] ++;

            if (i == n / 2 - 1) {
                for (int j = 2;j<=tot;j++) {
                    cnt2[j] = cnt[j];
                }
            }
        }

        for (int i = 2;i<=tot;i++) {
            cnt[i] -= cnt2[i];
        }
        return;
    }
};
using namespace PAM;

void solve() {
    int n;
    cin >> n >> s;
    s = s + s;  
    init();

    for (int i = 2;i<=tot;i++) {
        in[fail[i]] ++;
    }
    queue<int> q;
    for (int i = 2;i<=tot;i++) {
        if (in[i] == 0) q.push(i);
    }
    while (!q.empty()) {
        int u = q.front(); q.pop();
        cnt[fail[u]] += cnt[u];
        cnt[fail[u]] %= MOD;
        in[fail[u]] --;
        if (!in[fail[u]] && u > 1) {
            q.push(fail[u]);
        }
    }   

    ll ans = 0;
    for (int i = 2;i<=tot;i++) {
        if (len[i] > n) continue;
        ll res = 1ll * cnt[i] * cnt[i] % MOD;
        res = res * len[i] % MOD;
        ans = (ans + res) % MOD;
        // cout << cnt[i] << ' ' << len[i] << endl;    
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
