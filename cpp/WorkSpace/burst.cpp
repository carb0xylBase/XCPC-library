#include <bits/stdc++.h>
using namespace std;
#define DEBUG 1
#define endl '\n'
using ll = long long;
using pii = pair<int,int>;
using db = double;
using i128 = __int128_t;
const ll N = 2000000;
const db PI = acos(-1);
const db eps = 1e-9;

int c[N];

void init(int n) {
    for (int i = 1;i<=n;i++) c[i] = 0;
}

int lowbit(int x) { return x & -x; }

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

int main() {
#if DEBUG
    freopen("input.txt", "r", stdin);
    auto start_time = chrono::steady_clock::now();
#else
    ios::sync_with_stdio(false);
#endif
    cin.tie(nullptr);

    int n, k; if (!(cin >> n >> k)) return 0;
    vector<pair<int,int>> a(n);
    for (int i=0;i<n;i++) cin >> a[i].first >> a[i].second;
    vector<pair<int,int>> b(k);
    for (auto&v : b) {
        cin >> v.second >> v.first;
        v.second = PI - v.second;
        v.first = PI - v.first;
    }

    vector<int> ans(n, 0);
    for (auto [l, r] : b) {
        db sita = db(l + r) / 2 / 180 * PI;
        struct Node { db x,y; int idx; };
        vector<Node> nodes; nodes.reserve(n);
        for (int i = 0;i<n;i++) {
            auto [xx, yy] = a[i];
            Node z;
            z.x = xx * cos(sita) + yy * sin(sita);
            z.y = -xx * sin(sita) + yy * cos(sita);
            z.idx = i;
            nodes.push_back(z);
        }
        sita = db(r - l) / 2 / 180 * PI;
        db ksl = tan(sita);
        ksl = abs(ksl);
        auto update = [&](vector<db>& q) {
            sort(q.begin(), q.end());
            vector<db> nq;
            for (auto v : q) {
                if (nq.empty()) nq.push_back(v);
                else if (abs(nq.back() - v) > eps) nq.push_back(v);
            }
            swap(nq, q);
        };

        sort(nodes.begin(), nodes.end(), [&](const Node &L, const Node &R){
            db A = L.y + ksl * L.x;
            db B = R.y + ksl * R.x;
            if (abs(A - B) > eps) return A < B; // keyA ascending
            db a2 = L.y - ksl * L.x;
            db b2 = R.y - ksl * R.x;
            return a2 > b2; // keyB descending for tie on keyA
        });

        vector<db> q;
        q.reserve(n);
        for (auto &v : nodes) q.push_back(v.y - ksl * v.x);
        update(q);

        auto idx = [&](double x) -> int {
            int pos = lower_bound(q.begin(), q.end(), x) - q.begin();
            if (pos < (int)q.size() && abs(q[pos] - x) <= eps) return pos + 1;
            if (pos > 0 && abs(q[pos-1] - x) <= eps) return pos;
            return 0;
        };

        init(q.size());
        int total = 0;
        for (auto &v : nodes) {
            int id = idx(v.y - ksl * v.x);
            if (!id) continue;
            int less = (id - 1 >= 1) ? query(id - 1) : 0;
            ans[v.idx] += total - less;
            add(id, q.size());
            total++;
        }
    }

    for (auto v : ans) cout << v << ' ';
    cout << endl;

#if DEBUG
    auto end_time = chrono::steady_clock::now();
    auto diff = chrono::duration_cast<chrono::milliseconds>(end_time - start_time);
    cerr << "Time: " << diff.count() << " ms" << endl;
#endif
    return 0;
}
