#include <bits/stdc++.h>
using namespace std;
#define DEBUG 1
using ll = long long;
using pii = pair<int, int>;
using i128 = __int128_t;
using pll = pair<ll, ll>;
const ll N = 2000000;
const ll INF = 5e18;
const ll MOD = 1e9 + 7;

struct BigInt {
	vector<int> a;
	BigInt(int x = 0) {
		if (x == 0) {
			a.push_back(0);
			return;
		}
		while (x) {
			a.push_back(x % 10);
			x /= 10;
		}
	}
	void update() {
		for (int i = 0;i<a.size();i++) {
			if (a[i] >= 10) {
				if (i + 1 < a.size()) {
					a[i + 1] += a[i] / 10;
					a[i] %= 10;
				} else {
					a.push_back(a[i] / 10);
					a[i] %= 10;
				}
			}
		}
		while (a.size() > 1 && a.back() == 0) {
			a.pop_back();
		}
		return;
	}
	BigInt operator*(const BigInt& A) const {
		BigInt B;
		B.a.resize(a.size() + A.a.size());
		for (int i = 0;i<a.size();i++) {
			for (int j = 0;j<A.a.size();j++) {
				B.a[i + j] += a[i] * A.a[j];
			}
		}
		B.update();
		return B;
	}
};

ostream &operator<<(ostream &o, const BigInt &a) {
	for (int i = a.a.size() - 1;i>=0;i--) {
		o << a.a[i];
	}
	return o;
}