#include <stdio.h>
#include <cstdio>
#include <iostream>
#include <string.h>
#include <algorithm>
#include <queue>
#include <vector>
#include <stack>
#include <string>
#include <math.h>
#include <stdlib.h>
#include <deque>
#include <map>
#include <set>
#include <unordered_set>
#include <unordered_map>

#define M 1
#define MM 1
#define MMM 1
#define N 1
#define NN 4002

#define ll long long
#define ld double
#define vll vector<ll>
#define X first
#define Y second
#define PI 3.14159265358979323
#define FASTIO ios_base::sync_with_stdio(false); cin.tie(NULL); cout.tie(NULL);

#define fo(i,a,b) for(ll i = a; i <= b; i++)
#define foi(a) for(ll i=1;i<=a;i++)
#define foi0(a) for(ll i=0;i<a;i++)
#define foj(a) for(ll j=1;j<=a;j++)
#define foj0(a) for(ll j=0;j<a;j++)
#define fok(a) for(ll k=1;k<=a;k++)
#define fok0(a) for(ll k=0;k<a;k++)
#define fori for(ll i=1;i<=n;i++)
#define forin for(ll i=1;i<=n;i++)
#define forin1 for(ll i=n;i>=1;i--)
#define fori0 for(ll i=0;i<n;i++)
#define forj for(ll j=1;j<=m;j++)
#define forjn for(ll j=1;j<=n;j++)
#define forji for(ll j=1;j<=i;j++)
#define forjn0 for(ll j=0;j<n;j++)
#define forj0 for(ll j=0;j<m;j++)
#define fork for(ll k=1;k<=l;k++)
#define fork0 for(ll k=0;k<l;k++)
#define forkn for(ll k=1;k<=n;k++)
#define foriw for(ll i=1;;i++)

#define scans sc(s);slen=s.size()
#define scana fori sc(a[i]);
#define scanna scann; scana;
#define scanb fori sc(b[i])
#define scanbm forj sc(b[j])
#define scand fori sc(d[i])
#define scanaa fori forj sc(aa[i][j]);
#define scanaan fori forjn sc(aa[i][j]);
#define scanaa1 fori{scans;forj0{aa[i][j+1] = s[j] - '0';}s.clear();}
#define scanbb1 fori{scans;forj0{bb[i][j+1] = s[j] - '0';}s.clear();}
#define scanbb fori forj sc(bb[i][j]);
#define scanline(s) getline(cin,s); slen=s.size();
#define scanv v.resize(n); fori0 sc(v[i]);

#define prld(a,b) {cout << fixed; cout.precision(a); pr1(b);}
#define printcase pr("Case "); pr(++casenum);pr(": ")
#define printcases pr("Case #"); pr(++casenum);pr(": ")
#define printa fori {pr1(a[i]); }prl;
#define printal fori {pr1(a[i]); prl;}
#define printb fori {pr1(b[i]); }prl;
#define printd fori {pr1(d[i]); }prl;
#define printaa fori {forj {pr1(aa[i][j]);} prl;}prl;
#define printbb fori {forj {pr1(bb[i][j]);} prl;}prl;
#define printgg pr1l("gg");
#define printggg pr1l("ggg");
#define printv(v) for(auto qwe:v) {pr1(qwe);}; prl;
#define printvl(v) for(auto qwe:v) {pr1l(qwe);};
#define all(v) v.begin(), v.end()

#define sorta sort(a+1,a+n+1)
#define sortb sort(b+1,b+n+1)
#define sortd sort(d+1,d+n+1)
#define sortv(v) sort(all(v))
#define vsort(v) sort(v.begin(),v.end())
#define suma sum=0; fori sum+=a[i];

#define test pr1l("TEST!");
#define w1 while(1)
#define INF (ll) 1e18
#define boundcheck(tx,ty) if(tx>=1&&ty>=1&&tx<=n&&ty<=m)

#define X first
#define Y second
#define pb push_back

#define sc(a) cin >> a
#define sc1(a) cin >> a
#define sc2(a,b) cin >> a >> b
#define sc3(a,b,c) cin >> a >> b >> c
#define sc4(a,b,c,d) cin >> a >> b >> c >> d
#define sc5(a,b,c,d,e) cin >> a >> b >> c >> d >> e
#define sc6(a,b,c,d,e,f) cin >> a >> b >> c >> d >> e >> f

#define pr(a) cout << (a)
#define pr0 cout << (0)
#define prl cout << '\n'
#define pr1(a) cout << (a) << ' '
#define pr2(a,b) cout << (a) << ' ' << (b) << ' '
#define pr3(a,b,c) cout << (a) << ' ' << (b) << ' '<< (c) << ' '
#define pr4(a,b,c,d) cout << (a) << ' ' << (b) << ' '<< (c) << ' '<< (d) << ' '
#define pr5(a,b,c,d,e) cout << (a) << ' ' << (b) << ' '<< (c) << ' '<< (d) << ' '<< (e) << ' '
#define pr6(a,b,c,d,e,f) cout << (a) << ' ' << (b) << ' '<< (c) << ' '<< (d) << ' '<< (e) << ' ' << (f) << ' '
#define pr0l cout << '\n';
#define pr1l(a) cout << (a) << '\n'
#define pr2l(a,b) cout << (a) << ' ' << (b) << '\n'
#define pr3l(a,b,c) cout << (a) << ' ' << (b) << ' '<< (c) << '\n'
#define pr4l(a,b,c,d) cout << (a) << ' ' << (b) << ' '<< (c) << ' '<< (d) << '\n'
#define pr5l(a,b,c,d,e) cout << (a) << ' ' << (b) << ' '<< (c) << ' '<< (d) << ' '<< (e) << '\n'
#define pr6l(a,b,c,d,e,f) cout << (a) << ' ' << (b) << ' '<< (c) << ' '<< (d) << ' '<< (e) << ' ' << (f) << '\n'


using namespace std;
typedef pair<ll, ll> xy;
typedef vector<vll> matrix;
ll i, j, ii, jj, n, zz, yyy, xxx, maxim, ttttt, ja, mo, he, l1, l2, l3, mm, l4, end, zero, finish, tt, next, bre, cnt, ans, slen, to, casenum, nn, hab, count, t, now, one, two, yy, m, yes, cntt, x1, x2, x3, x4, y4, Y1, y2, y3, temp, i1, i2, J1, j2, i3, j3, len1, len2, low, mid, left, right, high, re, ok, last, tx, ty, k, num2, start, diff, cha, idx, num, xx, qq, w, e, no, r, sum, minim = INF, x, y, z, l, len, mini = INF, maxi = -INF;
ll dx[5] = { 0,0,1,0,-1 };
ll dy[5] = { 0,1,0,-1,0 };
ll ddx[9] = { 0,-1,-1,-1,0,0,1,1,1 };
ll ddy[9] = { 0,-1,0,1,-1,1,-1,0,1 };
ld ld1, ld0, ld2, ld3, ld4, ld5, ld6, ld7, lda[M], ldb[M];
ll a[N], d[N], b1[M], a1[M], a2[M], a3[M], a4[M], a5[M], bb[MM][MM];
ll b[N], tree[N], alis[M], dd[NN][NN], p[M], h[M], un[M], dist[M], aa[MM][MM], aa1[MM][MM], aa2[MM][MM], d1[M], d2[M];
ll dp[M][2], matn = 2, mu[M], tmp[N], suffix[N], aaa[MMM][MMM][MMM];
bool check[M], visit[N], treecheck[M], boo[M], visited[M], checkk[MM][MM];
char c1, c2, c, c3, c4, cc[M];
ld ldmax, ldmin, ldmax1, ldmax2, ldmin1, ldmin2, ldd[M];
ld ldx1,ldx2,ldx3,ldx4,ldy1,ldy2,ldy3,ldy4;

string str, s, s1, s2, s3, ss[M];
queue<ll> q, qx, qy;
priority_queue<ll, vector<ll>, greater<ll>> pq, pq2;
priority_queue<xy> pqxy;
priority_queue<xy> pqxy2;
stack<ll> st;
deque<ll> dq;
deque<xy> dqxy;
vll v[M], v1, v2, v3, print, scv[M], rscv[M], va[M];
vector<vll> vv, scc;
vector<xy> vxy, vxya[M], vxy2, vxy3;
vector<xy> vpa[M];
vector<pair<ll,xy>> vxyz;
xy xy1, xya[M];
ll mod = INF;

ll zegob(ll x, ll y){
    ll k = 1;
    while (y > 0) {
        if (y & 1) k = (k%mod*x%mod) % mod;
        x = (x*x) % mod;
        y >>= 1;
    }
    return k % mod;
}

bool vowel(char c){
    return (c == 'a' || c == 'e' || c == 'i' || c == 'o' || c == 'u');
}

bool bound(ll x, ll y, ll n, ll m){
    return (x > 0 && y > 0 && x <= n && y <= m);
}

ll gcd(ll x, ll y){
    if (x < y) swap(x, y);
    if (y == 0) return x;
    return gcd(y, x % y);
}

ll lcm(ll x, ll y){
    return x*y/gcd(x,y);
}

ll big(ll x, ll y){
    if (x > y) return x;
    return y;
}

ll small(ll x, ll y){
    if (x < y) return x;
    return y;
}

ld distance(ld x1, ld y1, ld x2, ld y2){
    return sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2));
}

ll lca(ll x, ll y) {

    if (d[x] > d[y])
        swap(x, y);
    for (i = 20; i >= 0; i--) {
        if (d[y] - d[x] >= (1 << i))
            y = dp[y][i];
    }
    if (x == y) return x;
    for (i = 20; i >= 0; i--) {
        if (dp[x][i] != dp[y][i]) {
            x = dp[x][i];
            y = dp[y][i];
        }
    }
    return dp[x][0];
}


void dfs(ll x, ll depth)
{
    visit[x] = true;
    d[x] = depth;
    ll l = v[x].size();
    for (ll i = 0; i < l;i++) {
        y = v[x][i];
        if (visit[y])
            continue;
        dp[y][0] = x;
        dfs(y, depth + 1);
    }
}

void solve() {
    maxi=0;
    sc(s1);
    sc(s2);
    n=s1.size();
    m=s2.size();
    fori0{
        forj0{
            if(s1[i]==s2[j]){
                dd[i+1][j+1]=dd[i][j]+1;
                maxi=max(maxi,dd[i+1][j+1]);
            }
        }
    };
    pr(maxi);
}

int main(void) {
    FASTIO;
    solve();
}
