#include <bits/stdc++.h>
using namespace std;
typedef long long ll;

const ll N = 2000000;
const ll MOD = 998244353;
const ll MAX = 1e18;

struct ShortestPath {
    int n;
    vector<vector<int>>& g;
    struct Node {
        ll val;
        int idx;
        bool operator<(const Node& A) const {

        };
    };
    priority_queue<Node> nodes;
    void init(vector<vector<int>>& G) {
        g = G;
        int n = g.size() + 1;
    }
    void cal() {
        while (!nodes.empty()) {
            Node u = nodes.top(); nodes.pop();
            for (auto v : g[u.idx]) {

            }
        }
    }
};

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
class Dijkstra{
public:
	Graph<T>* g;
	struct Node{
		int dist,idx;
		bool operator<(const Node &a) const {
			return dist > a.dist;
		}
		Node(int d,int x) {dist = d,idx = x;}
	};

	ll dist[N];
	bool vis[N];

	void dijkstra(int pos){
		priority_queue<Node> pq;
		for (int i = 1;i<g->n+1;i++){
			dist[i] = MAX;
			vis[i] = 0;
		}
		dist[pos] = 0;
		pq.push(Node(0,pos));

		while (!pq.empty()){
			int u = pq.top().idx;
			pq.pop();
			if (vis[u]){
				continue;
			}
			vis[u] = 1;
			for (int i = g->head[u];i;i = g->edge[i].next){
				int v = g->edge[i].to;
				if (dist[v] > dist[u] + g->edge[i].w){
					dist[v] = dist[u] + g->edge[i].w;
					if (!vis[v]){
						pq.push(Node(dist[v],v));
					}
				}
			}
		}
		return;
	}
};

/*
check 用于检查负环同时构造新边权
solve 输入起点,然后 dist 数组存最短路
*/

template<typename T>
class Johnson {
public:
    using ll = long long;
    Graph<T>* g;
    vector<T> h;
    vector<ll> dist;
    vector<bool> vis;

    void init(Graph<T>* G_) {
        g = G_;
        h.assign(g->n + 1, 0);
        vis.assign(g->n + 1, false);
    }

    bool check() {
        int n = g->n;
        int m = g->m;
        g->edge[0] = typename Graph<T>::Edge(); // dummy
        for (int i = 1; i <= n; i++) {
            g->add(0, i, 0); // 从虚拟源点 0 连向每个点
        }

        queue<int> q;
        vector<int> cnt(n + 2), inq(n + 2, 0);
        h.assign(n + 2, numeric_limits<T>::max() / 2);
        h[0] = 0;
        q.push(0);
        inq[0] = 1;

        while (!q.empty()) {
            int u = q.front(); q.pop();
            inq[u] = 0;
            for (int i = g->head[u]; i; i = g->edge[i].next) {
                int v = g->edge[i].to;
                T w = g->edge[i].w;
                if (h[v] > h[u] + w) {
                    h[v] = h[u] + w;
                    if (!inq[v]) {
                        q.push(v);
                        inq[v] = 1;
                        if (++cnt[v] > n) return false;
                    }
                }
            }
        }

        // 重新权重变换
        for (int u = 1; u <= n; u++) {
            for (int i = g->head[u]; i; i = g->edge[i].next) {
                int v = g->edge[i].to;
                g->edge[i].w += h[u] - h[v];
            }
        }

        return true;
    }

    void solve(int s) {
        dist.assign(g->n + 1, numeric_limits<ll>::max());
        vis.assign(g->n + 1, false);
        priority_queue<pair<ll, int>, vector<pair<ll, int>>, greater<>> q;

        dist[s] = 0;
        q.emplace(0, s);
        while (!q.empty()) {
            auto [d, u] = q.top(); q.pop();
            if (vis[u]) continue;
            vis[u] = true;
            for (int i = g->head[u]; i; i = g->edge[i].next) {
                int v = g->edge[i].to;
                T w = g->edge[i].w;
                if (dist[v] > dist[u] + w) {
                    dist[v] = dist[u] + w;
                    q.emplace(dist[v], v);
                }
            }
        }

        for (int i = 1; i <= g->n; i++) {
            if (dist[i] < numeric_limits<ll>::max() / 2) {
                dist[i] = dist[i] - h[s] + h[i];
            }
        }
    }
};