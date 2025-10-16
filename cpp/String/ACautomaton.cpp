#include <bits/stdc++.h>
using namespace std;
#define DEBUG 1
using ll = long long;
using pii = pair<ll,ll>;
const ll N = 2000000;
const ll INF = 5e18;
const ll MOD = 1e9 + 7;

namespace ACautomaton {
    constexpr int N = 2e6 + 6;
    constexpr int LEN = 2e6 + 6;
    constexpr int SIZE = 2e6 + 6;
    struct Node {
        int son[26],ans,fail,du,idx;

        void init() {
            memset(son,0,sizeof son);
            ans = fail = idx = 0;
        }
    } tr[SIZE];

    int tot,ans[N],pidx;

    void init() {
        tot = pidx = 0;
        tr[0].init();
    }

    void insert(string& s,int& idx) {
        int u = 0;
        for (int i = 0;i<s.size();i++) {
            int& son = tr[u].son[s[i]-'a'];
            if (!son) son = ++tot,tr[son].init();
            u = son;
        }

        if (!tr[u].idx) tr[u].idx = ++pidx;
        idx = tr[u].idx;
        return;
    }

    void build() {
        queue<int> q;
        for (int i = 0;i<26;i++) {
            if (tr[0].son[i]) q.push(tr[0].son[i]);
        }

        while (!q.empty()) {
            int u = q.front();
            q.pop();
            for (int i = 0;i<26;i++) {
                if (tr[u].son[i]) {
                    tr[tr[u].son[i]].fail = tr[tr[u].fail].son[i];
                    tr[tr[tr[u].fail].son[i]].du++;
                    q.push(tr[u].son[i]);
                } else {
                    tr[u].son[i] = tr[tr[u].fail].son[i];
                }
            }
        }
        return;
    }

    void query(string& t) {
        int u = 0;
        for (int i = 0;i<t.size();i++) {
            u = tr[u].son[t[i] - 'a'];
            tr[u].ans++;
        }
        return;
    } 

    void topu() {
        queue<int> q;
        for (int i = 0;i<=tot;i++) {
            if (tr[i].du == 0) q.push(i);
        }

        while (!q.empty()) {
            int u = q.front();
            q.pop();
            ans[tr[u].idx] = tr[u].ans;
            int v = tr[u].fail;
            tr[v].ans += tr[u].ans;
            if (!--tr[v].du) q.push(v);
        }

        return;
    }
};

int idx[N];
void solve() {
    int n; cin >> n;
    for (int i = 0;i<n;i++) {
        string s; cin >> s;
        ACautomaton::insert(s,idx[i]);
        ACautomaton::ans[i] = 0;
    }
    ACautomaton::build();
    string t; cin >> t;
    ACautomaton::query(t);
    ACautomaton::topu();

    int ANS = 0;
    for (int i = 0;i<n;i++) {
        if (ACautomaton::ans[idx[i]]) {
            ANS++;
        }
    }

    cout << ANS << endl;
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
