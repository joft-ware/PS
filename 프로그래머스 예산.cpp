#include <iostream>
#include <stdio.h>
#include <string>
#include <vector>
#include <algorithm>

using namespace std;

int solution(vector<int> d, int budget) {
    int answer = 0, sum=0;
    sort(d.begin(),d.end());
    for(auto i:d){
        if(i+sum>budget) break;
        sum+=i;
        answer++;
    }
    return answer;
}
