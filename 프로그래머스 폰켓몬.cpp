#include <vector>
#include <algorithm>
using namespace std;
bool mon[200001];

int solution(vector<int> nums)
{
    int answer = 0;
    for(auto i:nums){
        if(!mon[i]){
            mon[i]=true;
            answer++;
        }
    }
    if(answer>nums.size()/2)
        answer=nums.size()/2;
    return answer;
}
