#include <iostream>
#include <cstring>
using namespace std;

const int N=6000006,mod=998244353;
char s[N];

int pam[N][10],last,p,ni;
int fail[N],cnt[N],len[N],sz[N];
int lenlen;

int newnode(int l)
{
    memset(pam[p], 0, sizeof(pam[p]));
    cnt[p] = 0;
    // num[p] = 0;
    len[p] = l;
    return p ++;
}

void init()
{
    p = last = 0;
    newnode(0), newnode(-1);
    s[0] = '?';
    fail[0] = 1;
}

int get_fail(int x)
{
    while(s[ni - len[x] - 1] != s[ni])
        x = fail[x];
    return x;
}

void add(int c)
{
    int old = get_fail(last);
    if(!pam[old][c])
    {
        int now = newnode(len[old] + 2);
        fail[now] = pam[get_fail(fail[old])][c];
        pam[old][c] = now;
        // num[now] = num[fail[now]] + 1;
    }
    last = pam[old][c];
    cnt[last] ++;
}

void calc()
{
    for(int i = p-1; i>=0; i--)
    {
        cnt[fail[i]] += cnt[i];
    }
}

int main()
{
    freopen("input.txt", "r", stdin);
    init();
    
    scanf("%d", &lenlen);
    scanf("%s", s+1);
    for(int i=1+lenlen; i<=lenlen*2; i++)
        s[i] = s[i-lenlen];
    for(ni=1; ni<=lenlen; ni++)
        add(s[ni]-'0');
    for(int i=p-1; i>=0; i--)
        sz[i] = cnt[i];
    for(int i=p-1; i>=0; i--)
        sz[fail[i]] += sz[i];
    for(ni=lenlen+1; ni<=lenlen*2; ni++)
        add(s[ni] - '0');
    calc();
    long long ans = 0;
    for(int i=2; i<p; i++)
        if(len[i] <= lenlen)
            // cout << len[i] << " " << cnt[i] << "  " << sz[i] << "\n",
            ans = (ans + (1ll*(cnt[i] - sz[i])*(cnt[i] - sz[i])%mod) * len[i])%mod;
    cout << ans << '\n';
    return 0;
}