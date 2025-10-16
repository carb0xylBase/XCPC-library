#include <bits/stdc++.h>
using namespace std;
#define DEBUG 1
#define endl '\n'
#define FUCK if (DEBUG) cout << "fuck" << endl;
using ll = long long;
using pii = pair<int,int>;
using pll = pair<ll,ll>;
using db = double;
using pdd = pair<double, double>;
using i128 = __int128_t;
const ll N = 2000000;
const ll INF = 1e6;
const ll MOD = 1e9 + 7;
const db PI = acos(-1);
const db eps = 1e-9;

int c[N];

void init(int n) {
    for (int i = 1;i<n+1;i++) c[i] = 0;
}

int lowbit(int x) {
    return x & -x;
}

void add(int x, int n) {
    while (x <= n) {
        c[x] ++;
        x += lowbit(x);
    }
}

int query(int x) {
    int ans = 0;
    while (x) {
        ans += c[x];
        x -= lowbit(x);
    }
    return ans;
}
void solve() {
    int n, k; cin >> n >> k;
    vector<pii> a(n);
    for (auto&v : a) cin >> v.first >> v.second;
    vector<pii> seg(k);
    for (auto&v : seg) cin >> v.first >> v.second;
    vector<int> ans(n, 0);

    struct BIT {
        int n; vector<int> t;
        void init(int _n){ n=_n; t.assign(n+1,0); }
        void add(int i,int v=1){ for(; i<=n; i+=i&-i) t[i]+=v; }
        int query(int i){ int r=0; for(; i>0; i-=i&-i) r+=t[i]; return r; }
    };

    for (auto [L, R] : seg) {
        double mid = double(L + R) / 2.0 / 180.0 * PI;
        double half = double(R - L) / 2.0 / 180.0 * PI;
        double kf = fabs(tan(half));

        struct Node {
            double u; // rotated x
            double q1; // k*u - v
            double p2; // v + k*u
            int idx;
        };
        vector<Node> nodes; nodes.reserve(n);
        for (int i = 0; i < n; ++i) {
            double x = a[i].first, y = a[i].second;
            double ux = x * cos(mid) + y * sin(mid);
            double uy = -x * sin(mid) + y * cos(mid);
            nodes.push_back({ux, kf * ux - uy, uy + kf * ux, i});
        }

        // compress p2
        vector<double> p2s; p2s.reserve(n);
        for (auto &nd : nodes) p2s.push_back(nd.p2);
        sort(p2s.begin(), p2s.end());
        p2s.erase(unique(p2s.begin(), p2s.end(), [&](double A, double B){ return fabs(A-B) <= eps; }), p2s.end());
        auto get_p2_id = [&](double x)->int{
            return int(upper_bound(p2s.begin(), p2s.end(), x + eps) - p2s.begin());
        };

        // sort by u ascending, and if equal, sort by q1 then p2 to make grouping contiguous
        sort(nodes.begin(), nodes.end(), [&](const Node &A, const Node &B){
            if (fabs(A.u - B.u) > eps) return A.u < B.u;
            if (fabs(A.q1 - B.q1) > eps) return A.q1 < B.q1;
            return A.p2 < B.p2;
        });

        // CDQ on array [0..n-1], we ensure we never split equal u across halves
        vector<Node> tmp(n);
        BIT bit; bit.init((int)p2s.size()+5);

        function<void(int,int)> cdq = [&](int l, int r){
            if (l >= r) return;
            int mid = (l + r) >> 1;
            // move mid right to avoid splitting equal u's
            while (mid < r && fabs(nodes[mid].u - nodes[mid+1].u) <= eps) ++mid;
            if (mid == r) { // all equal u in [l,r], no pairs with u_j < u_i
                return;
            }
            cdq(l, mid);
            cdq(mid+1, r);

            // merge by q1: left part [l..mid], right part [mid+1..r]
            int i = l, j = mid+1, t = l;
            // we will add left nodes into BIT when left.q1 <= right.q1
            vector<int> added_ids;
            while (j <= r) {
                while (i <= mid && nodes[i].q1 <= nodes[j].q1 + eps) {
                    int id = get_p2_id(nodes[i].p2);
                    bit.add(id, 1);
                    added_ids.push_back(id);
                    ++i;
                }
                int idr = get_p2_id(nodes[j].p2);
                ans[nodes[j].idx] += bit.query(idr);
                ++j;
            }
            // clear BIT additions
            for (int id : added_ids) bit.add(id, -1);

            // now standard merge to keep nodes sorted by q1 for upper levels (not strictly necessary but keep for stability)
            i = l; j = mid+1; t = l;
            while (i <= mid && j <= r) {
                if (nodes[i].q1 <= nodes[j].q1 + eps) tmp[t++] = nodes[i++]; else tmp[t++] = nodes[j++];
            }
            while (i <= mid) tmp[t++] = nodes[i++];
            while (j <= r) tmp[t++] = nodes[j++];
            for (int p = l; p <= r; ++p) nodes[p] = tmp[p];
        };

        cdq(0, n-1);
    }

    for (int i = 0; i < n; ++i) {
        cout << ans[i] << (i+1==n?'\n':' ');
    }
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