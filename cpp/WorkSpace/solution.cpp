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
const ll MOD = 1e9 + 7;

class Graph {
public:
    int n;        
    struct Edge {
        int next, to, rev;
        ll w;
    } edge[2 * N];  
    int head[N];      
    int cnt; 

    void init(int N_) {
        n = N_;
        for (int i = 0; i <= n; i++)
            head[i] = 0;
        cnt = 0;
    }

    void add_edge(int u, int v, ll w) {
        cnt++;
        edge[cnt].to = v;
        edge[cnt].w = w;
        edge[cnt].next = head[u];
        head[u] = cnt;
        int forward_index = cnt;
        cnt++;
        edge[cnt].to = u;
        edge[cnt].w = 0;
        edge[cnt].next = head[v];
        head[v] = cnt;
        int reverse_index = cnt;
        edge[forward_index].rev = reverse_index;
        edge[reverse_index].rev = forward_index;
    }
} g;

class Dinic {
private:
    Graph* G;       
    int s, t;             
    vector<ll> d;         
    vector<int> cur;       

    bool bfs() {
        fill(d.begin(), d.end(), -1);
        queue<int> q;
        d[s] = 0;
        q.push(s);
        while (!q.empty()) {
            int u = q.front(); q.pop();
            for (int i = G->head[u]; i; i = G->edge[i].next) {
                int v = G->edge[i].to;
                if (d[v] < 0 && G->edge[i].w > 0) {
                    d[v] = d[u] + 1;
                    q.push(v);
                }
            }
        }
        return d[t] >= 0;
    }

    ll dfs(int u, ll flow) {
        if (u == t) return flow;
        for (int &i = cur[u]; i; i = G->edge[i].next) {
            int v = G->edge[i].to;
            if (d[v] == d[u] + 1 && G->edge[i].w > 0) {
                ll pushed = dfs(v, min(flow, G->edge[i].w));
                if (pushed) {
                    G->edge[i].w -= pushed;
                    G->edge[G->edge[i].rev].w += pushed;
                    return pushed;
                }
            }
        }
        return 0;
    }

public:
    Dinic(Graph* graph, int source, int sink) : G(graph), s(source), t(sink) {
        d.resize(G->n + 1);
        cur.resize(G->n + 1);
    }

    ll max_flow() {
        ll flow = 0;
        while (bfs()) {
            for (int i = 0; i <= G->n; i++)
                cur[i] = G->head[i];
            while (ll pushed = dfs(s, INF))
                flow += pushed;
        }
        return flow;
    }
};

void solve() {
    int n; cin >> n;
    vector<int> a(n + 1, 0), b(n + 1, 0);
    for (int i = 1;i<n+1;i++) cin >> a[i];
    for (int i = 1;i<n+1;i++) cin >> b[i];
    vector<bool> vis(n + 1, 0);
    vector<vector<int>> gp[2];
    for (int i = 1;i<n+1;i++) {
        if (vis[i]) continue;
        vector<int> q;
        int u = i;
        while (!vis[u]) {
            vis[u] = 1;
            u = a[u];
            q.push_back(u);
        }
        gp[0].push_back(q);
    }
    for (int i = 1;i<n+1;i++) {
        vis[i] = 0;
    }
    for (int i = 1;i<n+1;i++) {
        if (vis[i]) continue;
        vector<int> q;
        int u = i;
        while (!vis[u]) {
            vis[u] = 1;
            u = b[u];
            q.push_back(u);
        }
        gp[1].push_back(q);
    }

    // for (auto v : gp[1]) {
    //     for (auto w : v) {
    //         cout << w << ' ';
    //     }
    //     cout << endl;
    // }
    // cout << endl;

    int S = gp[0].size() + gp[1].size() + 2 * n + 1;
    int E = S + 1;
    g.init(E);
    for (int i = 0;i<gp[0].size();i++) {
        g.add_edge(S, 2 * n + i + 1, 1);
        for (auto v : gp[0][i]) {
            g.add_edge(2 * n + i + 1, v, 1);
        }
    }
    for (int i = 0;i<gp[1].size();i++) {
        g.add_edge(2 * n + gp[0].size() + i + 1, E, 1);
        for (auto v : gp[1][i]) {
            g.add_edge(v + n, 2 * n + gp[0].size() + i + 1, 1);
        }
    }
    for (int i = 1;i<=n;i++) {
        g.add_edge(i, i + n, 1);
    }

    Dinic solver(&g, S, E);
    int z = solver.max_flow();
    vector<int> ans;
    for (int i = 1;i<=n;i++) {
        for (int j = g.head[i];j;j=g.edge[j].next) {
            int v = g.edge[j].to;
            int w = g.edge[j].w;
            if (v == i + n && w == 1) {
                ans.push_back(i);
                break;
            }
        }
    }
    cout << ans.size() << endl;
    for (auto v : ans) {
        cout << v << ' ';
    }
    cout << endl;
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
