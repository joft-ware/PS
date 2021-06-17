#include <string>
#include <vector>
#include <algorithm>

using namespace std;

vector<int> solution(vector<int> array, vector<vector<int>> commands) {
    vector<int> answer;
    int i, j, k;
    for(auto command:commands){
        vector<int> temp=array;
        i=command[0]-1;
        j=command[1]-1;
        k=command[2]-1;
        sort(temp.begin()+i,temp.begin()+j+1);
        answer.push_back(temp[i+k]);
    }
    return answer;
}
