#include <string>
#include <vector>
#include <algorithm>
#include <stdio.h>

using namespace std;

vector<int> solution(vector<int> answers) {
    vector<int> answer;
    int index=0, maxi;
    int a1[5]={5,1,2,3,4};
    int a2[8]={5,2,1,2,3,2,4,2};
    int a3[10]={5,3,3,1,1,2,2,4,4,5};
    int arr[4]={0};
    for(auto i: answers){
        index++;
        if(a1[index%5]==i) arr[1]++;
        if(a2[index%8]==i) arr[2]++;
        if(a3[index%10]==i) arr[3]++;
    }
    maxi = max(arr[1],max(arr[2],arr[3]));    
    for(int i=1;i<=3;i++) if(arr[i]==maxi) answer.push_back(i);    
    return answer;
}

