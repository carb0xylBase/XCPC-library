#include <bits/stdc++.h>
using namespace std;
#define DEBUG 1
#define FUCK if (DEBUG) cout << "fuck" << endl;
using ll = long long;
using pii = pair<int,int>;
using i128 = __int128_t;
using pll = pair<ll,ll>;
const ll N = 2000000;
const ll INF = 1e18;
const ll MOD = 1e9 + 7;

namespace ST {
    const int M = 21;
    int lg[N];

    void init(int n) {
        lg[0] = 0, lg[1] = 0;
        for (int i = 2;i<=n;i++) {
            lg[i] = lg[i / 2] + 1;
        }
    }

    struct ST {
        struct Info {
            ll mx;

            Info (ll MX = 0) {
                mx = MX;
            }

            Info operator+(const Info& A) {
                return Info(max(mx, A.mx));
            }
        };

        vector<vector<Info>> f;

        void init(vector<ll>& a) {
            int n = a.size() - 1;
            f.clear();
            f.resize(n + 1);
            for (int i = 1;i<=n;i++) {
                f[i].resize(M + 1);
            }
            for (int j = 1;j<=M;j++) {
                for (int i = 1;i+(1<<j)-1<=n;i++) {
                    f[i][j] = f[i][j - 1] + f[i + (1 << (j - 1))][j - 1];
                }
            }
        }

        Info query(int l, int r) {
            if (l > r) return Info();
            int s = lg[r - l + 1];
            return f[l][s] + f[r - (1 << s) + 1][s];
        }
    };
};
