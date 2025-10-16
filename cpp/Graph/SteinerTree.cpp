#include <bits/stdc++.h>
using namespace std;

typedef long long ll;
const ll N = 2000000;

template<typename T>
class Graph{
public:
    int n,m;
    struct Edge{
        int next,to;
        T w;
    };
    vector<Edge> edge;
    vector<int> head;
    int cnt;
    void add(int u,int v,T w){
        edge[++cnt].next = head[u];
        head[u] = cnt;
        edge[cnt].to = v;
        edge[cnt].w = w;
    }
	void init(int N,int M){
		n = N; m = M;
		head.clear();
        head.resize(n+1,0);
        Edge.clear();
        Edge.resize(M + 1);
		cnt = 0;
		return;
	}
};

template<typename T>
class SteinerTree {
public:
    Graph<T>* g;
    int n, k;
    vector<int> spe;
    vector<vector<T>> dp;
    vector<T> dist;
    vector<bool> vis;

    void init(Graph<T>* _g, const vector<int>& _spe) {
        g = _g; n = _g->n; k = _spe.size(); spe = _spe;
        int U = 1 << k;
        dp.assign(n+1, vector<T>(U, numeric_limits<T>::max()/4));
        dist.assign(n+1, numeric_limits<T>::max()/4);
        vis.assign(n+1, false);
        for (int i = 0; i < k; i++)
            dp[spe[i]][1<<i] = 0;
    }

    T solve() {
        int U = (1<<k) - 1;
        for (int S = 1; S <= U; S++) {
            for (int A = (S-1)&S; A; A = (A-1)&S)
                for (int i = 1; i <= n; i++)
                    dp[i][S] = min(dp[i][S], dp[i][A] + dp[i][S^A]);

            priority_queue<pair<T,int>, vector<pair<T,int>>, greater<>> pq;
            for (int i = 1; i <= n; i++) {
                dist[i] = dp[i][S];
                vis[i] = false;
                pq.emplace(dist[i], i);
            }
            while (!pq.empty()) {
                auto [d,u] = pq.top(); pq.pop();
                if (vis[u]) continue;
                vis[u] = true;
                for (int e = g->head[u]; e; e = g->edge[e].next) {
                    int v = g->edge[e].to; T w = g->edge[e].w;
                    if (dist[v] > d + w) {
                        dist[v] = d + w;
                        pq.emplace(dist[v], v);
                    }
                }
            }
            for (int i = 1; i <= n; i++)
                dp[i][S] = dist[i];
        }

        T ans = numeric_limits<T>::max()/4;
        for (int i = 1; i <= n; i++)
            ans = min(ans, dp[i][U]);
        return ans;
    }
};