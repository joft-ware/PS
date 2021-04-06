#include <stdio.h>
#include <cstdio>
#include <iostream>
#include <algorithm>
#include <queue>
#include <vector>
using namespace std;
priority_queue<int, vector<int>, greater<int>> pq;
int i, n, t;
int main(void) {
    ios_base::sync_with_stdio(false); cin.tie(NULL); cout.tie(NULL);
    cin >> n;
    for(i=1;i<=n*n;i++){
        cin >> t;
        if(pq.size()<n) pq.push(t);
        else if(pq.top()<t){
            pq.pop();
            pq.push(t);
        }
    }
    cout << (pq.top());
}