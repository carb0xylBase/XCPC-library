#include <bits/stdc++.h>
using namespace std;

using ll = long long;

class NumberTheoryBlock {
public:
    void down(ll x,vector<pair<ll,ll>>& block) {
        block.clear();
        ll l = 1,r = 0;
        while (l <= x) {
            r = x / (x / l);
            if (r > x) r = x;
            block.emplace_back(l,r);
            l = r + 1;
        }
        return;
    }

    void up(ll x,vector<pair<ll,ll>>& block) {
        block.clear();
        ll l = 1,r = 0;
        while (l <= x) {
            if (l == x) {
                block.emplace_back(l,l);
                break;
            }
            r = (x - 1) / ((x - 1) / l);
            block.emplace_back(l,r);
            l = r + 1;
        }
        return;
    }
};