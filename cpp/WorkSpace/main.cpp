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
const ll INF = 1e18;
const ll MOD = 1e9 + 7;

namespace ST {
    int lg[N];
    void init(int n) {
        lg[0] = 0, lg[1] = 0;
        for (int i = 2;i<=n;i++) {
            lg[i] = lg[i / 2] + 1;
        }
    }
    struct ST {
        struct Info {
            ll val;
            Info() {
                val = INF;
            }
            Info operator+(const Info& A) const {
                Info z;
                z.val = min(val, A.val);
                return z;
            }
        };
        vector<vector<Info>> f;
        void init(const vector<ll>& a) {
            int n = (int)a.size();
            if (n == 0) return;
            for (int i = 2; i <= n; ++i) lg[i] = lg[i >> 1] + 1;
            int K = lg[n];
            f.assign(K + 1, vector<Info>(n));
            for (int i = 0; i < n; ++i) f[0][i].val = a[i];
            for (int k = 1; k <= K; ++k)
                for (int i = 0; i + (1 << k) <= n; ++i)
                    f[k][i] = f[k - 1][i] + f[k - 1][i + (1 << (k - 1))];
        }
        Info query(int l, int r) {
            if (l > r) return Info();
            int len = r - l + 1;
            int k = lg[len];
            return f[k][l] + f[k][r - (1 << k) + 1];
        }
    } pp, ss;
}


struct SegTree {
    #define LS (rt << 1)
    #define RS (rt << 1 | 1)
    struct Node {
        ll mn;
        Node() {
            mn = INF;
        }
    } nodes[N * 4];

    Node merge(Node L, Node R) {
        Node z;
        z.mn = min(L.mn, R.mn);
        return z;
    }

    void build(int rt, int l, int r) {
        if (l == r) {
            nodes[rt].mn = INF;
            return;
        }
        int mid = l + r >> 1;
        build(LS, l, mid);
        build(RS, mid + 1, r);
        nodes[rt] = merge(nodes[LS], nodes[RS]);
    }

    Node query(int rt, int l, int r, int ql, int qr) {
        if (r < ql || qr < l || ql > qr) {
            return Node();
        }
        if (ql <= l && r <= qr) {
            return nodes[rt];
        }
        int mid = l + r >> 1;
        return merge(query(LS, l, mid, ql, qr), query(RS, mid + 1, r, ql, qr));
    }

    void modfiy(int rt, int l, int r, int q, ll val) {
        if (q < l || q > r || l > r) return;
        if (l == r) {
            nodes[rt].mn = val;
            return;
        }
        int mid = l + r >> 1;
        if (q <= mid) {
            modfiy(LS, l, mid, q, val);
        } else {
            modfiy(RS, mid + 1, r, q, val);
        }
        nodes[rt] = merge(nodes[LS], nodes[RS]);
        return;
    }
} pre, suf;

void solve() {
    int n, k; cin >> n >> k;
    vector<ll> a(n + 1, 0);
    for (int i = 1;i<n+1;i++) cin >> a[i];
    string s; cin >> s;

    pre.modfiy(1, 0, n + 1, 0, 0);
    suf.modfiy(1, 0, n + 1, n + 1, 0);
    int lst = 0;
    for (int i = 1;i<=n;i++) {
        ll val = pre.query(1, 0, n + 1, max(lst, i - k), i - 1).mn;
        val += a[i];
        pre.modfiy(1, 0, n + 1, i, val);
        if (s[i - 1] == '1') {
            lst = i;
        }
    }
    lst = n + 1;
    for (int i = n;i>=1;i--) {
        ll val = suf.query(1, 0, n + 1, i + 1, min(lst, i + k)).mn;
        val += a[i];
        suf.modfiy(1, 0, n + 1, i, val);
        if (s[i - 1] == '1') {
            lst = i;
        }
    }

    vector<ll> P, S;
    for (int i = 0;i<=n+1;i++) {
        P.push_back(pre.query(1, 0, n + 1, i, i).mn);
    }
    for (int i = 0;i<=n+1;i++) {
        S.push_back(suf.query(1, 0, n + 1, i, i).mn);
    }
    ST::pp.init(P);
    ST::ss.init(S);

    int q; cin >> q;
    while (q--) {
        ll p, v; cin >> p >> v;
        ll ans = INF;
        ll L = max(0ll, p - k), R = min(p + k, 1ll * n + 1);
        for (int i = p - 1;i>=L;i--) {
            if (s[i - 1] == '1') {
                L = i;
                break;
            }
        }
        for (int i = p + 1;i<=R;i++) {
            if (s[i - 1] == '1') {
                R = i;
                break;
            }
        }
        ans = ST::pp.query(L, p - 1).val
            + ST::ss.query(p + 1, R).val + v;
        // ans = pre.query(1, 0, n + 1, L, p - 1).mn + 
        //     suf.query(1, 0, n + 1, p + 1, R).mn + v;
        
        if (s[p - 1] != '1') {
            for (int i = L;i<p;i++) {
                int r = min({R, n + 1ll, i + k * 1ll});
                ans = min(ans, ST::pp.query(i, i).val + 
                    ST::ss.query(p + 1, r).val);
                // ans = min(ans, pre.query(1, 0, n + 1, i, i).mn + 
                //     suf.query(1, 0, n + 1, p + 1, r).mn);
            }
        }

        cout << ans << endl;
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
    ST::init(6e5);
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
