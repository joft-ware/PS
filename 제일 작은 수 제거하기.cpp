#include <string>
#include <vector>
#include <algorithm>

using namespace std;

vector<int> solution(vector<int> arr) {
    vector<int> answer;
    int min=arr[0];
    for(auto i:arr) if(min>i) min=i;
    for(auto i:arr) if(min!=i) answer.push_back(i);
    if(arr.size()==1) answer.push_back(-1);
    return answer;
}
