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

#define M 100001
#define MM 101
#define N 100001

#define ll long long
#define ull unsigned ll
#define ld double
#define vll vector<ll>
#define vs vector<string>
#define Yes "Yes"
#define No "No"
#define YES "YES"
#define NO "NO"
#define mn m=n;
#define X first
#define Y second
#define PI 3.14159265358979323846264338327950288419716939937510

#define foi(a) for(ll i=1;i<=a;i++)
#define foi0(a) for(ll i=0;i<a;i++)
#define foj(a) for(ll j=1;j<=a;j++)
#define foj0(a) for(ll j=0;j<a;j++)
#define fok(a) for(ll k=1;k<=a;k++)
#define fori for(ll i=1;i<=n;i++)
#define forin for(ll i=1;i<=n;i++)
#define forin1 for(ll i=n;i>=1;i--)
#define fori0 for(ll i=0;i<n;i++)
#define forj for(ll j=1;j<=m;j++)
#define forjn for(ll j=1;j<=n;j++)
#define forji for(ll j=1;j<i;j++)
#define forjn0 for(ll j=0;j<n;j++)
#define forj0 for(ll j=0;j<m;j++)
#define fork for(ll k=1;k<=l;k++)
#define fork0 for(ll k=0;k<l;k++)
#define forkn for(ll k=1;k<=n;k++)
#define foriw for(ll i=1;;i++)

#define scann sc(n)
#define scanm sc(m);
#define scant sc(t);
#define scanx sc(x);
#define scany sc(y);
#define scank sc(k);
#define scanc sc(c);
#define scanxy sc2(x,y);
#define scanyx scanf("%lld %lld",&y,&x)
#define scanzr scanf("%lld %lld",&z,&r)
#define scanwe scanf("%lld %lld",&w,&e)
#define scannm sc2(n,m);
#define scanwe scanf("%lld %lld",&w,&e)
#define scanmn scanf("%lld %lld",&m,&n)
#define scannml scanf("%lld %lld %lld",&n,&m,&l)
#define scanxyz sc3(x,y,z);
#define scanxyzr scanf("%lld %lld %lld %lld",&x,&y,&z,&r)
#define scans sc(s);slen=s.size();
#define scansn frees; scanf("%s", &s[1]); len = strlen(&s[1]); n=len;
#define scansm frees; scanf("%s", &s[1]); len = strlen(&s[1]); m=len;
#define scans1 cin >> s1; len1 = s1.size();
#define scans2 cin >> s1; len1 = s1.size();
#define scana freea; fori sc(a[i]);
#define scanna scann; fori scanf("%lld",&a[i]);
#define scana1d fori scanf("%1d",&a[i]);
#define scanb1d fori scanf("%1d",&b[i]);
#define scanb fori scanf("%lld",&b[i]);
#define scand fori scanf("%lld",&d[i]);
#define scanaa fori for(ll j=1;j<=m;j++) scanf("%lld",&aa[i][j]);
#define scanaa1 fori{scans;forj0{aa[i][j+1] = s[j] - '0';}s.clear();}
#define scanbb1 fori{scans;forj0{bb[i][j+1] = s[j] - '0';}s.clear();}
#define scanbb fori for(ll j=1;j<=m;j++) scanf("%lld",&bb[i][j]);
#define scanline(s) getline(cin,s); slen=s.size();
#define scanv vll v(n+1); fori sc(v[i]);
#define scannv scann; scanv;

#define prld(a,b) {cout << fixed; cout.precision(a); pr1(b);}
#define printld(a) prld(a)
#define printsum printf("%lld\n",sum);
#define printcase pr("Case "); pr(++casenum);pr(": ");
#define printcases pr("Case #"); pr(++casenum);pr(": ");
#define prints printf("%s",&s[1]);
#define printc printf("%c",c);
#define printmax printf("%lld\n",maxi);
#define printmin printf("%lld\n",mini);
#define printmini printf("%lld\n",mini);
#define printa fori {printf("%lld ",a[i]); }printf("\n");
#define printa1 fori {printf("%lld ",a1[i]); }printf("\n");
#define printa2 fori {printf("%lld ",a2[i]); }printf("\n");
#define printb fori {printf("%lld ",b[i]); }printf("\n");
#define printd fori {printf("%lld ",d[i]); }printf("\n");
#define printaa fori {for(ll j=1;j<=m;j++) {printf("%3d ",aa[i][j]);} printf("\n");}printf("\n");
#define printbb fori {for(ll j=1;j<=m;j++) {printf("%3d ",bb[i][j]);} printf("\n");}printf("\n");
#define printgg pr1l("gg");
#define printv(v) for(auto qwe:v) {pr1(qwe);}; prl;
#define prv(v) printv(v)
#define printvl(v) for(auto qwe:v) {pr1l(qwe);};
#define prvl(v) printvl(v)
#define pra printa
#define prcase printcase
#define prcases printcases

#define frees for(ll i=0;i<=len+n;i++) s[i]=0;
#define freea for(ll i=0;i<=n;i++) a[i]=0;
#define cleana for(ll i=0;i<=n;i++) a[i]=0;
#define cleanb for(ll i=0;i<=n;i++) b[i]=0;
#define cleanaa(x) for(ll iii=0;iii<=n+1;iii++) {for(ll jjj=0;jjj<=m+1;jjj++) {aa[iii][jjj]=x;}};
#define sorta sort(a+1,a+n+1);
#define sortb sort(b+1,b+n+1);
#define sortd sort(d+1,d+n+1);
#define sortv(v) sort(full(v));
#define suma sum=0; fori sum+=a[i];
#define infa fori a[i]=INF;
#define reversea fori tempa[i]=a[n+1-i]; fori a[i]=tempa[i];
#define findmax maxi=a[1]; fori if(a[i]>maxi) maxi=a[i];
#define findmaxn maxi=a[1]; fori if(a[i]>maxi) {maxi=a[i]; num=i;}
#define findmin mini=a[1]; fori if(a[i]<mini) mini=a[i];
#define issmall(a)      ((a>='a')&&(a<='z'))
#define isbig(a)      ((a>='A')&&(a<='Z'))

#define lens len = strlen(s);
#define test printf("TEST!");
#define wt while(t--)
#define w1 while(1)
#define INF (ll)0x7f7f7f7f
#define br break
#define braek break
#define bk break
#define nbreak if(n==0) break;
#define boundcheck(tx,ty) if(tx>=1&&ty>=1&&tx<=n&&ty<=m)
#define strint fori a[i] = s[i] - '0';
#define full(v) v.begin(), v.end()
#define all(v) v.begin(), v.end()

#define X first
#define Y second
#define fo(i,a,b) for(ll i = a; i <= b; i++)
#define pb push_back
#define mp make_pair
#define pbm(a,b) push_back(make_pair(a,b))
#define vsort(v) sort(v.begin(),v.end());
#define sc(a) cin >> a
#define sc1(a) cin >> a
#define sc2(a,b) cin >> a >> b
#define sc3(a,b,c) cin >> a >> b >> c
#define sc4(a,b,c,d) cin >> a >> b >> c >> d
#define sc5(a,b,c,d,e) cin >> a >> b >> c >> d >> e
#define sc6(a,b,c,d,e,f) cin >> a >> b >> c >> d >> e >> f
#define scn sc(n)
#define scm sc(m)
#define scnm sc2(n,m)
#define scx sc(x)
#define sct sc(t)
#define scxy sc2(x,y)
#define scline scanline

#define pr(a) cout << (a)
#define pr0 cout << (0);
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
#define pr7l(a,b,c,d,e,f,g) cout << (a) << ' ' << (b) << ' '<< (c) << ' '<< (d) << ' '<< (e) << ' ' << (f) << ' ' << (g) << '\n'
#define pr8l(a,b,c,d,e,f,g,h) cout << (a) << ' ' << (b) << ' '<< (c) << ' '<< (d) << ' '<< (e) << ' ' << (f) << ' ' << (g) << ' ' << (h) << '\n'

#define prcnt pr1l(cnt)
#define prx pr1l(x)
#define prxy pr2l(x, y)
#define prno pr1l("no")
#define pryes pr1l("yes")
#define prNO pr1l("NO")
#define prYES pr1l("YES")
#define prgg pr1l("gg")
#define prmaxi pr1l(maxi)
#define prmax pr1l(maxi)
#define prmin pr1l(mini)
#define prmini pr1l(mini)
#define prnum pr1l(num)
#define prsum printsum
#define prstr for(ll wq=1;wq<=slen;wq++) pr(str[wq]);

using namespace std;
ll i, j, ii, jj, n, zz, yyy, xxx, maxim, ja, mo, l1, l2, l3, mm, l4, end, zero, finish, next, bre, cnt, ans, slen, to, casenum, nn, hab, count, t, now, one, two, yy, m, yes, cntt, x1, x2, x3, x4, y4, Y1, y2, y3, temp, i1, i2, J1, j2, i3, j3, len1, len2, low, mid, left, right, high, ok, tx, ty, k, num2, start, num, xx, qq, w, e, no, r, sum, minim = INF, x, y, z, l, len, mini = INF, maxi = -INF, x11, x22, x33, y11, y22, y33;
ll dx[5] = { 0,0,1,0,-1 };
ll dy[5] = { 0,1,0,-1,0 };
ll ddx[9] = { 0,-1,-1,-1,0,0,1,1,1 };
ll ddy[9] = { 0,-1,0,1,-1,1,-1,0,1 };
ll dddx[9] = { 0,-1,-2,-1,0,0,1,2,1 };
ll dddy[9] = { 0,-1,0,1,-2,2,-1,0,1 };
ll knightdx[9] = { 0,-1,-1,1,1,-2,-2,2,2 };
ll knightdy[9] = { 0,2,-2,2,-2,1,-1,-1,1 };
ll alphabet_lines[27] = {0,3,2,1,2,4,3,1,3,1,1,3,1,3,2,1,2,2,2,1,2,1,1,1,2,2,1};
ld ld1, ld2, ld3, ld4, ld5, ld6, ld7, lda[M], ldb[M];
ll a[M], b1[M], a1[M], a2[M], a3[M], a4[M], a5[M], rank[M], bb[MM][MM], sumtree[M], mintree[M], maxtree[M], minindextree[M], prime[M];
ll b[M], alis[M], dd[MM][MM], p[M], h[M], ax[M], ay[M], az[M], d[M], dist[M], aa[MM][MM], d1[M], d2[M], tempa[M], lazy[M];
ll qry[M][4], dp[MM][11], matn = 2, mu[M];
bool check[M], visit[M], treecheck[M], boo[M];
char c1, c2, c, c3, c4, cc[M];
ld ldmax, ldmin, ldmax1, ldmax2, ldmin1, ldmin2, ldd[M];

string str, s, s1, s2, s3, ss[M], ss1[M], ss2[M];
typedef pair<ll, ll> xy;
typedef vector<vll> matrix;
ull u1, u2, u3, u4;
queue<ll> q, qx, qy;
priority_queue<ll> pq, pq2;
priority_queue<xy> pqxy;
stack<ll> st;
deque<ll> dq;
deque<xy> dqxy;
vll v, v1, v2, v3, print, vv[M], rv[M], visited;
vector<xy> vxy, vxya[M];
xy xy1, xya[M];
vector<vll> scc;
ll mod = INF;

matrix operator *(matrix &a, matrix &b) {
    matrix c(2, vll(2));
    ll n = m = l = 2;
    fori0
        forj0
            fork0
                c[i][j] = (c[i][j] + (a[i][k] * b[k][j])) % mod;
    return c;
}

ll zegobmod(ll x, ll y)
{
    ll k = 1;
    while (y > 0) {
        if (y & 1)
            k = (k%mod*x%mod) % mod;
        x = (x*x) % mod;
        y >>= 1;
    }
    return k % mod;
}

ll zegob(ll x, ll y)
{
    ll k = 1;
    while (y > 0) {
        if (y & 1)
            k = (k*x);
        x = (x*x);
        y >>= 1;
    }
    return k;
}

ll binary_search(vll v, ll x, ll left, ll right) {
    ll mid = (left + right) / 2;
    if (v[mid] == x)
        return mid;
    if (left >= right)
        return -1;
    if (v[mid] < x)
        return binary_search(v, x, mid + 1, right);
    return binary_search(v, x, left, mid);
}

bool da(char c)
{
    if (c >= 'A' && c <= 'Z')
        return true;
    return false;
}

bool so(char c)
{
    if (c >= 'a' && c <= 'z')
        return true;
    return false;
}

char to_lower(char c){
    if(c<'a')
        return c+32;
    return c;
}

char to_upper(char c){
    if(c>='a')
        return c-32;
    return c;
}

bool isnum(char c) {
    return (c >= '0'&&c <= '9');
}

bool daso(char c) {
    return (da(c) || so(c));
}

bool vowel(char c)
{
    if (c == 'a' || c == 'e' || c == 'i' || c == 'o' || c == 'u')
        return true;
    return false;
}

bool bound(ll x, ll y, ll n, ll m)
{
    if (x > 0 && y > 0 && x <= n && y <= m)
        return true;
    return false;
}

ll find(ll x) {// p:부모
    if (p[x] == x)
        return x;
    p[x] = find(p[x]);
    return p[x];
}

void uni(ll x, ll y) {// h: height
    x = find(x);
    y = find(y);
    if (x == y) return;
    if (h[x] > h[y])
        swap(x, y);
    p[x] = y;
    if (h[x] == h[y]) h[y]++;
}

bool same(ll x, ll y)
{
    return (find(x) == find(y));
}

ll gcd(ll x, ll y)
{
    if (x < y)
        swap(x, y);
    if (y == 0)
        return x;
    return gcd(y, x % y);
}

ll lcm(ll x, ll y) {
    return x * y / gcd(x, y);
}

ll bigger(ll x, ll y)
{
    if (x > y)
        return x;
    return y;
}

ll smaller(ll x, ll y)
{
    if (x < y)
        return x;
    return y;
}

ll find_max(long long* a, ll n)
{
    findmax;
    return maxi;
}

void clean(long long* a, int n)
{
    fori
        a[i] = 0;
    a[0] = 0;
}

ll zari(ll n)
{
    k = 10;
    foriw
    {
        if (k > n)
            return i;
        k *= 10;
    }
}

ll biggest(ll x, ll y, ll z)
{
    ll a[4];
    a[1] = x;
    a[2] = y;
    a[3] = z;
    sort(a + 1, a + 4);
    return a[3];
}

ll smallest(ll x, ll y, ll z)
{
    ll a[4];
    a[1] = x;
    a[2] = y;
    a[3] = z;
    sort(a + 1, a + 4);
    return a[1];
}

ll minindex(ll x, ll y) {
    if (a[x] == a[y]) return (smaller(x, y));
    return (a[x] < a[y]) ? x : y;
}

ll maketree_minindex(ll left, ll right, ll node)
{
    if (left >= right)
        return minindextree[node] = left;
    else
    {
        ll mid = (left + right) / 2;
        ll leftnode = maketree_minindex(left, mid, node * 2);
        ll rightnode = maketree_minindex(mid + 1, right, node * 2 + 1);
        return minindextree[node] = minindex(leftnode, rightnode);
    }
}

ll query_minindex(ll node, ll left, ll right, ll start, ll end)
{
    if (right < start || end < left)
        return 0; // 겹치지 않는 경우(영향이 없는 값을 반환)
    if (start <= left && right <= end) return minindextree[node]; // 모두 겹치는 경우
    int mid = (left + right) / 2; // 일부만 겹치는 경우
    ll leftnode = query_minindex(node * 2, left, mid, start, end);
    ll rightnode = query_minindex(node * 2 + 1, mid + 1, right, start, end);
    return minindex(leftnode, rightnode);
}

ll maketree_min(ll left, ll right, ll node)
{
    if (left == right)
        return mintree[node] = a[left];
    else
    {
        int mid = (left + right) / 2;
        mintree[node] = smaller(maketree_min(left, mid, node * 2), maketree_min(mid + 1, right, node * 2 + 1)); //작은거
        return mintree[node];
    }
}

ll query_min(ll node, ll left, ll right, ll start, ll end)
{
    if (right < start || end < left)
        return INF; // 겹치지 않는 경우(영향이 없는 값을 반환)
    if (start <= left && right <= end) return mintree[node]; // 모두 겹치는 경우
    int mid = (left + right) / 2; // 일부만 겹치는 경우
    return smaller(query_min(node * 2, left, mid, start, end), query_min(node * 2 + 1, mid + 1, right, start, end));
}

ll update_min(ll node, ll left, ll right, ll idx) { ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    if (idx<left || idx>right) // update 필요없음
        return mintree[node];
    if (left == right)
        return mintree[node] = a[left];

    ll mid = (left + right) / 2;
    ll leftnode = update_min(node * 2, left, mid, idx);
    ll rightnode = update_min(node * 2 + 1, mid + 1, right, idx);
    mintree[node] = min(leftnode, rightnode);
    return mintree[node];
}

ll maketree_max(ll left, ll right, ll node)
{
    if (left == right)
        return maxtree[node] = a[left];
    else
    {
        ll mid = (left + right) / 2;
        maxtree[node] = bigger(maketree_max(left, mid, node * 2), maketree_max(mid + 1, right, node * 2 + 1)); //작은거
        return maxtree[node];
    }
}

ll query_max(ll node, ll left, ll right, ll start, ll end)
{
    if (right < start || end < left)
        return -1; // 겹치지 않는 경우(영향이 없는 값을 반환)
    if (start <= left && right <= end) return maxtree[node]; // 모두 겹치는 경우
    int mid = (left + right) / 2; // 일부만 겹치는 경우
    return bigger(query_max(node * 2, left, mid, start, end), query_max(node * 2 + 1, mid + 1, right, start, end));
}

ll maketree_sum(ll left, ll right, ll node)
{
    if (left == right)
    {
        sumtree[node] = a[left];
        return sumtree[node];
    }
    else
    {
        ll mid = (left + right) / 2;
        sumtree[node] = (maketree_sum(left, mid, node * 2) + maketree_sum(mid + 1, right, node * 2 + 1));
        return sumtree[node];
    }
}

void update_lazy_sum(ll node, ll left, ll right) {
    if (!lazy[node])
        return;
    sumtree[node] += ((right - left + 1)*lazy[node]);
    if (right != left) {
        lazy[node * 2] += lazy[node];
        lazy[node * 2 + 1] += lazy[node];
    }
    lazy[node] = 0;
}

ll update_sum(ll left, ll right, ll val, ll node, ll start, ll end) {
    update_lazy_sum(node, left, right);
    if (end<left || start>right) // 범위 밖
        return sumtree[node];
    if (start <= left && end >= right) { // 범위 내부에 속함
        lazy[node] += val;
        update_lazy_sum(node, left, right);
        return sumtree[node];
    }
    ll mid = (left + right) / 2; // 걸쳐 있음
    sumtree[node] = update_sum(left, mid, val, node * 2, start, end) + update_sum(mid + 1, right, val, node * 2 + 1, start, end);
    return sumtree[node];
}

ll query_sum(ll node, ll left, ll right, ll start, ll end) {
    update_lazy_sum(node, left, right);
    if (right < start || end < left)
        return 0; // 겹치지 않는 경우(영향이 없는 값을 반환)
    if (start <= left && end >= right) // 범위 내부에 속함
        return sumtree[node]; // 포함되는 경우
    ll mid = (left + right) / 2; // 일부만 겹치는 경우
    return (query_sum(node * 2, left, mid, start, end)) + (query_sum(node * 2 + 1, mid + 1, right, start, end));
}

ll fact(ll n)
{
    ll k = 1;
    fori
        k *= i;
    return k;
}

ll reverse(ll x)
{
    ll sum = 0;
    ll t = x;
    ll z = zari(x);
    foi(z)
    {
        sum += (t % 10) * zegob(10, z - i);
        t /= 10;
    }
    return sum;
}

ll ds(char c)
{
    ll x;
    if (so(c))
        x = (c - 'a' + 1);
    else
        x = (c - 'A' + 27);
    return x;
}

ll ab(ll x)
{
    if (x < 0)
        return -x;
    return x;
}

ll find_sum(ll left, ll right, ll node, ll sum) {
    update_lazy_sum(node, left, right);
    if (left >= right)
        return right;
    mid = (left + right) / 2;
    if (sumtree[node * 2] >= sum)
        return find_sum(left, mid, node * 2, sum);
    else
        return find_sum(mid + 1, right, node * 2 + 1, sum - sumtree[node * 2]);
}

ll insert_sum(ll node, ll left, ll right, ll start, ll end) {
    update_lazy_sum(node, left, right);

    mid = (left + right) / 2;
    if (left == right)
        return left;
    if (sumtree[node * 2] > 0 && mid >= start)
        return insert_sum(node * 2, left, mid, start, end);
    else
        return insert_sum(node * 2 + 1, mid + 1, right, start, end);
}

vll getpi(string p) { // 문자열 p의 pi배열 반환
    ll n = (ll)p.size();
    ll j = 0;
    vll pi(n, 0);
    fo(i, 1, n - 1) {
        while (j > 0 && p[i] != p[j]) // i: 기준, j: 비교 대상
            j = pi[j - 1];
        if (p[i] == p[j])
            pi[i] = ++j;
    };
    return pi;
}

vll kmp(string s, string s2) { // 문자열 s에 문자열 s2가 포함된 위치 벡터를 반환
    vll ans;
    auto pi = getpi(s2);
    ll n = (ll)s.size(), m = (ll)s2.size(), j = 0;
    fori0{
        while (j > 0 && s[i] != s2[j])
            j = pi[j - 1];
        if (s[i] == s2[j]) {
            if (j == m - 1) {
                ans.pb(i - m + 1);
                j = pi[j];
            }
            else
                j++;
        }
    };
    return ans;
}

vll getpi(vll p) { // 벡터 p의 pi배열 반환
    ll n = (ll)p.size();
    ll j = 0;
    vll pi(n, 0);
    fo(i, 1, n - 1) {
        while (j > 0 && p[i] != p[j]) // i: 기준, j: 비교 대상
            j = pi[j - 1];
        if (p[i] == p[j])
            pi[i] = ++j;
    };
    return pi;
}

vll kmpll(vll v1, vll v2) { // 벡터 v1에 벡터 v2가 포함된 위치 벡터를 반환
    vll ans;
    auto pi = getpi(v2);
    ll n = (ll)v1.size(), m = (ll)v2.size(), j = 0;
    fori0{
        while (j > 0 && v1[i] != v2[j])
            j = pi[j - 1];
        if (v1[i] == v2[j]) {
            if (j == m - 1) {
                ans.pb(i - m + 1);
                j = pi[j];
            }
            else
                j++;
        }
    };
    return ans;
}

ll lis(ll a[], ll n) { // lis 길이 구하기
    vll v;
    ll cnt = 0;
    fori alis[i] = 0;
    v.pb(a[1]);
    fo(i, 2, n) {
        if (v[cnt] < a[i])
        {
            v.pb(a[i]);
            alis[i] = ++cnt;
        }
        else {
            x = lower_bound(v.begin(), v.end(), a[i]) - v.begin();
            v[x] = a[i];
            alis[i] = x;
        }
    }
    return cnt + 1;
}

vll lisv(ll a[], ll n) { // a의 LIS 벡터 구하기
    vll v, ret;
    ll cnt = 0;
    fori d[i] = 0;
    fori visit[i] = false;
    v.pb(a[1]);
    fo(i, 2, n) {
        if (v[cnt] < a[i])
        {
            v.pb(a[i]);
            d[i] = ++cnt;
        }
        else {
            x = lower_bound(full(v), a[i]) - v.begin();
            v[x] = a[i];
            d[i] = x;
        }
    }
    visit[cnt + 1] = true;
    for (i = n; i >= 1; i--)
    {
        if (visit[d[i] + 1] && !visit[d[i]])
        {
            visit[d[i]] = true;
            ret.pb(i);
        }
    }
    return ret;
}

vll ntov(ll n) { // 정수 n을 vector 로 변환
    vll a;
    w1{
        a.pb(n % 10);
        if (n < 10)
            break;
        n /= 10;
    }
    return a;
}

vll stov(string s) //문자열 s를 vector로 변환
{
    vll a;
    ll n=s.size();
    fori0 a.pb(s[i]-'0');
    return a;
}

ll banolim(ld a) {
    ll x = (ll)a;
    ld y = (ld)x;
    if (a - y >= 0.5)
        return x + 1;
    else
        return x;
}

vll dijk(vector<xy> vpa[], ll start, ll n) { // 다익스트라. vpa: {to, cost}
    fori d[i] = INF;
    fori check[i] = false;
    priority_queue<xy> ppq;
    vll v;
    d[start] = 0;

    ppq.push({ -d[start],start }); // cost, 위치
    while (!ppq.empty()) {
        ll now = ppq.top().second;
        ll cost = -(ppq.top().first);
        ppq.pop();

        if (check[now]) continue;
        check[now] = true;

        ll l = vpa[now].size();
        fo(i, 0, l - 1) {
            ll y = vpa[now][i].first;
            ll tcost = vpa[now][i].second;
            if (check[y]) continue;
            if (d[y] > cost + tcost)
            {
                d[y] = cost + tcost;
                ppq.push({ -d[y],y });
            }
        };
    }

    fori v.pb(d[i]);
    return v;
}

ld ccw(ld x1, ld x2, ld x3, ld y1, ld y2, ld y3) {
    ld x = (x1*y2 + x2 * y3 + x3 * y1);
    x += (-y1 * x2 - y2 * x3 - y3 * x1);
    return x / 2;
}

ld ccw(xy a, xy b, xy c) {
    ld w = (b.X - a.X)*(c.Y - a.Y) - (b.Y - a.Y)*(c.X - a.X);
    if (w < 0) return -1;
    return (w > 0);
}

bool cross(xy a, xy b, xy c, xy d) { // 선분ab와 cd의 cross 여부
    ll x = ccw(a, b, c)*ccw(a, b, d);
    ll y = ccw(c, d, a)*ccw(c, d, b);
    if (!x && !y) {
        if (a > b) swap(a, b);
        if (c > d) swap(c, d);
        if (a <= d && b >= c) return true;
        return false;
    }
    return (x <= 0 && y <= 0);  // 등호 붙이면 접하는 경우 포함. 안붙이면 불포함
}

ld distxy(xy a, xy b) { // 좌표 거리
    ld w = a.X - b.X;
    ld e = a.Y - b.Y;
    return sqrt(w*w + e * e);
}

bool ccwcmp(xy a, xy b) {
    if (ccw(xy1, a, b) < 0) return true;
    if (ccw(xy1, a, b) > 0) return false;
    if (distxy(xy1, a) < distxy(xy1, b)) return true;
    return false;
}

bool xycmp(xy a, xy b) { // 시계방향
    if (a.Y > b.Y)
        return true;
    if (a.Y < b.Y)
        return false;
    return (a.X < b.X); // >
}

vector<xy> convex_hull(xy xya[], ll n) {
    vector<xy> vxy;
    xy1 = xya[1];
    fori if (xycmp(xya[i], xy1)) xy1 = xya[i]; // 극값 검색
    sort(xya + 1, xya + n + 1, ccwcmp);
    foi(n) { // 시계방향
        while (vxy.size() >= 2 && ccw(vxy[vxy.size() - 2], vxy[vxy.size() - 1], xya[i]) >= 0)
        {
            vxy.pop_back();
        }
        vxy.pb(xya[i]);
    };
    return vxy;
}

bool xycmpmax(xy a, xy b) {
    if (a.Y > b.Y)
        return true;
    if (a.Y < b.Y)
        return false;
    return (a.X > b.X);
}

bool xycmpmin(xy a, xy b) {
    if (a.Y < b.Y)
        return true;
    if (a.Y > b.Y)
        return false;
    return (a.X < b.X);
}

ld rotating_calipers(vector<xy> vxy) { // 시계방향으로 회전하는 캘리퍼스
    xy mini = vxy[0];
    xy maxi = vxy[0];
    ll x = 0, y = 0;
    ll n = vxy.size();
    fori0{
        if (xycmpmin(vxy[i],mini)) {
            mini = vxy[i];
            x = i;
        }
        if (xycmpmax(vxy[i],mini)) {
            maxi = vxy[i];
            y = i;
        }
    };
    l1 = vxy[x].X;
    l2 = vxy[x].Y;
    l3 = vxy[y].X;
    l4 = vxy[y].Y;
    ld maxim = distxy(vxy[x], vxy[y]);
    foi0(n) {
        ll nextx = (x + 1) % n;
        ll nexty = (y + 1) % n;
        xy xx = { vxy[x].X - vxy[nextx].X,vxy[x].Y - vxy[nextx].Y };
        xy yy = { vxy[y].X - vxy[nexty].X,vxy[y].Y - vxy[nexty].Y };
        if (ccw(xx, { 0,0 }, yy) < 0) x = nextx; // 원점으로 두 벡터를 옮긴 뒤 비교
        else y = nexty;
        if (maxim < distxy(vxy[x], vxy[y])) {
            maxim = max(maxim, distxy(vxy[x], vxy[y]));
            l1 = vxy[x].X;
            l2 = vxy[x].Y;
            l3 = vxy[y].X;
            l4 = vxy[y].Y;
        }
    };
    return maxim;
}

vll changebase(ll n, ll m) {
    vll a;
    vll b;
    w1{
        a.pb(n%m);
        cnt++;
        if (n < m)
            break;
        n /= m;
    };
    foi0(cnt) {
        if (a[cnt-1 - i] <= 9)
            b.pb(a[cnt-1 - i] + '0');
        else
            b.pb(a[cnt-1 - i] + 55);
    }
    return b;
}

ll fibo(ll n) {
    matrix x = { {1, 0},
                 {0, 1} };
    matrix a = { {1, 1},
                 {1, 0} };
    while (n > 0) {
        if (n % 2)
            x = x * a;
        a = a * a;
        n /= 2;
    }
    return x[0][1] % mod;
}

ll fibosum(ll from, ll to) {
    ll x = fibo(from + 1) - 1;
    ll y = fibo(to + 2) - 1;
    return (y - x + mod) % mod;
}

void binary_search(void) {
    l = 0;
    r = a[n];
    w1{
        mid = (l + r) / 2;
        sum = 1;
        x = a[1];
        fori{
            if (a[i] - x >= mid)
            {
                sum++;
                x = a[i];
            }
        }
        if (sum < m) { // 왼쪽(값을 작게 해야 함)
            if (l == r) break;
            r = mid;
        }
        else { // 오른쪽(값을 크게 해야 함)
            maxi = max(maxi,mid);
            if (l == r) break;
            l = mid + 1;
        }
        if (l > r) break;
    };
    pr1(maxi);
}

void union_find(void) { // 최적화된 union find (균형 트리)
    fori h[i] = 1, p[i] = i; // h:높이, p: 부모
}

struct student {
    string name;
    ll kor, eng, math;
};

bool cmp(student a, student b) {
    if (a.kor != b.kor) return (a.kor < b.kor);
    if (a.eng != b.eng) return(a.eng > b.eng);
    if (a.math != b.math) return (a.math < b.math);
    return a.name > b.name;
}

string lcs(string st1, string st2) {
    ll maxi = 0, cnt = 0;
    ll dd[MM][MM] = { { 0 } };
    string s3;
    stack<ll> st;
    s1 = '0' + st1;
    s2 = '0' + st2;
    ll n = s1.size();
    ll m = s2.size();
    fori0{
        forj0{
            if (i == 0 || j == 0) { dd[i][j] = 0; continue; }
            if (s1[i] == s2[j]) dd[i][j] = dd[i - 1][j - 1] + 1;
            else dd[i][j] = bigger(dd[i - 1][j],dd[i][j - 1]);
        }
    }
    ll i = n - 1;
    ll j = m - 1;
    while (dd[i][j] != 0) {
        if (dd[i][j] == dd[i][j - 1])
            j--;
        else if (dd[i][j] == dd[i - 1][j])
            i--;
        else {
            st.push(i);
            i--;
            j--;
        }
    }
    while (!st.empty()) {
        s3 += s1[st.top()];
        st.pop();
    }
    return s3;
}

void scc_(ll x, vll &list) {
    visit[x]=true;
    for(auto y:vv[x]) if(!visit[y]) scc_(y, list);
    list.pb(x);
}

void scc_reverse(ll x, vll &list) {
    visit[x]=true;
    for(auto y:rv[x]) if(!visit[y]) scc_reverse(y, list);
    list.pb(x);
}

void make_mu(){
    mu[0]=0;
    mu[1]=1;
    foi(N){
        for(ll j=i*2;j<=N;j+=i)
            mu[j]-=mu[i];
    }
}

ll square_free(ll x){
    ll cnt = 0;
    for(ll i=1;i*i<=x;i++) cnt+=mu[i]*x/(i*i);
    return cnt;
}

int main(void) {
    // 속도 최적화 // --------------------------------------------------
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);
    cout.tie(NULL);
    // 속도 최적화 // --------------------------------------------------
    scn;
    auto v=changebase(n,2);
    foi0(v.size()) sum+=zegob(3,v.size()-1-i)*(v[i]-'0');
    prsum;
};