#include <bits/stdc++.h>
using namespace std;
typedef long long ll;
const ll N = 2000000;
const ll INF = 1e18;

#define LS (rt << 1)
#define RS (rt << 1 | 1)

class SegTree{
public:
    struct Node{
        Node() {

        }
    }nodes[N*4];

    Node merge(Node L,Node R){
        Node M;
        return M;
    }
    void build(int rt,int l,int r){
        if (l == r){
            return;
        }
        int mid = l + r >> 1;
        build(LS,l,mid),build(RS,mid+1,r);
        nodes[rt] = merge(nodes[LS],nodes[RS]);
    }
    void pd(int rt){

    }
    void update(int rt,int l,int r,int ql,int qr,int val){
        if (ql > qr || l > qr || ql > r) return;
        if (ql <= l && r <= qr){
            return;
        }
        int mid = l+r>>1;
        pd(rt);
        if (ql <= mid){
            update(LS,l,mid,ql,qr,val);
        }
        if (qr >= mid + 1){
            update(RS,mid+1,r,ql,qr,val);
        }
        nodes[rt] = merge(nodes[LS],nodes[RS]);
        return;
    }
    void modify(int rt,int l,int r,int q,int val){
        if (l > q || r < q) return;
        if (l == r) {
            return;
        }
        int mid = l+r>>1;
        pd(rt);
        if (q <= mid){
            modify(LS,l,mid,q,val);
        }
        if (q >= mid + 1){
            modify(RS,mid+1,r,q,val);
        }
        nodes[rt] = merge(nodes[LS],nodes[RS]);
        return;
    }
    Node query(int rt,int l,int r,int ql,int qr){
        if (ql > qr) {
            return Node();
        }
        if (ql <= l && r <= qr){
            return nodes[rt];
        }
        int mid = l+r>>1;
        pd(rt);
        if (ql > mid){
            return query(RS,mid+1,r,ql,qr);
        }else if (qr < mid + 1){
            return query(LS,l,mid,ql,qr);
        }else{
            return merge(query(LS,l,mid,ql,qr),query(RS,mid+1,r,ql,qr));
        }
    }
};

class SegTreePlusMul{
struct Node{
    ll sum,plus,mul,len;
};
private:
    Node nodes[N];
public:
    ll mod;
    Node merge(Node A,Node B){
        Node C;
        C.plus = 0;
        C.mul = 1;
        C.sum = (A.sum + B.sum) % mod;
        C.len = A.len + B.len;
        return C;
    }
    void build(ll rt,ll l,ll r){
        if (l==r){
            cin >> nodes[rt].sum;
            nodes[rt].len = 1;
            nodes[rt].plus = 0;
            nodes[rt].mul = 1;
            return;
        }
        ll mid = l+r>>1;
        build(LS,l,mid),build(RS,mid+1,r);
        nodes[rt] = merge(nodes[LS],nodes[RS]);
    }
    void pushDown(int rt){
        nodes[RS].sum = (nodes[RS].sum*nodes[rt].mul%mod+nodes[rt].plus*nodes[RS].len%mod) % mod;
        nodes[LS].sum = (nodes[RS].sum*nodes[rt].mul%mod+nodes[rt].plus*nodes[LS].len%mod) % mod;

        nodes[RS].mul = (nodes[RS].mul*nodes[rt].mul) % mod;
        nodes[LS].mul = (nodes[LS].mul*nodes[rt].mul) % mod;

        nodes[RS].plus = (nodes[RS].plus*nodes[rt].mul+nodes[rt].plus) % mod;
        nodes[LS].plus = (nodes[LS].plus*nodes[rt].mul+nodes[rt].plus) % mod;

        nodes[rt].plus = 0;
        nodes[rt].mul = 1;
    }
    void update(ll rt,ll l,ll r,ll ql,ll qr,ll mode,ll k){
        if (ql <= l && r <= qr){
            if (mode == 1){
                nodes[rt].plus += k;
                nodes[rt].plus %= mod;
                nodes[rt].sum += (k * nodes[rt].len) % mod;
                nodes[rt].sum %= mod;
            }else{
                nodes[rt].plus *= k;
                nodes[rt].plus %= mod;
                nodes[rt].mul *= k;
                nodes[rt].mul %= mod;
                nodes[rt].sum *= k;
                nodes[rt].sum %= mod;
            }
            return;
        }
        pushDown(rt);
        ll mid = l+r>>1;
        if (ql<=mid){
            update(LS,l,mid,ql,qr,mode,k);
        }
        if (qr >= mid+1){
            update(RS,mid+1,r,ql,qr,mode,k);
        }
        nodes[rt] = merge(nodes[LS],nodes[RS]);
    }
    Node query(ll rt,ll l,ll r,ll ql,ll qr){
        if (ql <= l && r <= qr){
            return nodes[rt];
        }
        pushDown(rt);
        ll mid = l+r>>1;
        if (qr<=mid){
            return query(LS,l,mid,ql,qr);
        }else if (ql>=mid+1){
            return query(RS,mid+1,r,ql,qr);
        }else{
            return merge(query(LS,l,mid,ql,qr),query(RS,mid+1,r,ql,qr));
        }
    }
};


// 调用前先 init
#define LS nodes[rt].ls
#define RS nodes[rt].rs
struct DynamicSegTree {
    struct Node {
        int ls,rs;
        Node() {
            ls = rs = 0;
        }
    } nodes[N];
    int cntNode;
    void init() {
        cntNode = 1;
        nodes[0] = Node();
    }
    int newNode() {
        nodes[cntNode] = Node();
        return cntNode++;
    }
    Node merge(Node L,Node R) {
        Node M;
        // 补充合并方法
        return M;
    }
    void pushUp(int rt) {
        // 更新父节点
        return;
    }
    void pushDown(int rt) {
        // 下传标记
    }
    void update(int& rt,int l,int r,int ql,int qr,int val) {
        if (ql > qr || l > qr || r < ql) return;
        if (!rt) rt = newNode();
        if (ql <= l && r <= qr) {
            // 更新方法
            return;
        }
        int mid = l + r >> 1;
        update(nodes[rt].ls,l,mid,ql,qr,val);
        update(nodes[rt].rs,mid+1,r,ql,qr,val);
        pushUp(rt);
        return;
    }
    Node query(int& rt,int l,int r,int ql,int qr) {
        if (ql > qr || l > qr || r < ql || !rt) return Node();
        if (ql <= l && r <= qr) {
            return nodes[rt];
        }
        int mid = l + r >> 1;
        pushDown(rt);
        return merge(query(LS,l,mid,ql,qr),query(RS,mid+1,r,ql,qr));
    }
    int treeMerge(int& ls,int& rs,int l,int r) {
        if (!ls || !rs) return ls | rs;
        if (l == r) {
            // 叶节点合并方法
            
            return ls;
        }
        int mid = l + r >> 1;
        nodes[ls].ls = treeMerge(nodes[ls].ls,nodes[rs].ls,l,mid);
        nodes[ls].rs = treeMerge(nodes[ls].rs,nodes[rs].rs,mid+1,r);
        pushUp(ls);
        return ls;
    }
};


struct DynamicSegTree {
	#define LS info[rt].ls
	#define RS info[rt].rs
	struct Info {
        int ls, rs;
        Info() {
            ls = rs = 0;
        }
        Info operator+(const Info& A) const {

        }
    } info[N];
	int cntNode = 0;
	
	void init() {
		cntNode = 0;
	};

	int newNode() {
		++cntNode;
		info[cntNode] = Info();
		return cntNode;
	};

	void pushUp(int rt) {
		return;
	}

	void build(int& rt,int l,int r) {
		rt = newNode();
		if (l == r) {
			return;
		}
		int mid = l + r >> 1;
		build(LS,l,mid); build(RS,mid+1,r);
		pushUp(rt);
		return;
	}

    // 主席树版本
	// void update(int& rt,int l,int r,int ql,int qr) {
	// 	if (ql > qr || qr < l || r < ql) return;
	// 	if (ql <= l && r <= qr) {
	// 		return;
	// 	}
	// 	int mid = l + r >> 1;
	// 	if (ql <= mid) {
	// 		++cntNode;
	// 		info[cntNode] = info[LS];
	// 		LS = cntNode;
	// 		update(LS,l,mid,ql,qr);
	// 	}
	// 	if (qr >= mid + 1) {
	// 		++cntNode;
	// 		info[cntNode] = info[RS];
	// 		RS = cntNode;
	// 		update(RS,mid+1,r,ql,qr);
	// 	}
	// 	pushUp(rt);
	// 	return;
	// }

    void update(int& rt,int l,int r,int ql,int qr) {
		if (ql > qr || qr < l || r < ql) return;
        if (!rt) rt = newNode();
		if (ql <= l && r <= qr) {
			return;
		}
		int mid = l + r >> 1;
		if (ql <= mid) {
			update(LS,l,mid,ql,qr);
		}
		if (qr >= mid + 1) {
			update(RS,mid+1,r,ql,qr);
		}
		pushUp(rt);
		return;
	}

	Info query(int rt,int l,int r,int ql,int qr) {
		if (!rt || ql > qr || r < ql || qr < l) return Info();
		if (ql <= l && r <= qr) {
			return info[rt];
		}
		int mid = l + r >> 1;
		if (qr <= mid) {
			return query(LS,l,mid,ql,qr);
		} else if (ql >= mid + 1) {
			return query(RS,mid+1,r,ql,qr);
		} else {
			return query(LS,l,mid,ql,qr) + query(RS,mid+1,r,ql,qr);
		}
	}
} solver;

