#include <vector>
#include <iostream>
using namespace std;
int a[3001];
int solution(vector<int> nums) {
    int answer = 0;
    int i, j, k, len=nums.size();
    for(int i=2;i<=3000;i++){
        if(a[i]==1) continue;
        a[i]=2;
        for(int j=i*2;j<=3000;j+=i) a[j]=1;
    }        
    for(i=0;i<len;i++)
        for(j=i+1;j<len;j++)
            for(k=j+1;k<len;k++)
                if(a[nums[i]+nums[j]+nums[k]]==2) answer++;
    return answer;
}
