#define _USE_MATH_DEFINES
#include <bits/stdc++.h>
using namespace std;
typedef unsigned long long ULL;

#include "../derived_once_nets.cpp"


int findNewId(const vector<MarkedCombinatorialMap>& old_list,
              const vector<MarkedCombinatorialMap>& new_list,
              int old_id) {
    assert(0 <= old_id and old_id < ssize(old_list));
    const MarkedCombinatorialMap& A = old_list[old_id];
    auto B = symmetrize(A);
    for (int new_id = 0; new_id < ssize(new_list); new_id++) {
        if (areIsomorphic(A, new_list[new_id])
            or areIsomorphic(B, new_list[new_id])) {
            return new_id;
        }
    }
    return -36178914;
}

int main() {
    const vector<MarkedCombinatorialMap>& old_list = parseMarkedCombinatorialMaps("../data/old_derived_once_computer.txt");
    const vector<MarkedCombinatorialMap>& new_list = parseMarkedCombinatorialMaps("../data/derived_once_computer.txt");

    vector<int> old_ids = {
        28, 100, 160, 184, 252, 322, 328, 350, 365, 426, 474, 490, 526, 621, 690, 722, 747, 790, 839, 924, 935, 954, 963, 993, 1048, 1118, 1122, 1147, 1160,100,101,11,39,722, 360,45,1122,1067};

    cout << "old -> new" << endl;
    for (int old_id: old_ids) {
        int new_id = findNewId(old_list, new_list, old_id);
        cout << old_id << " -> " << new_id << endl;
    }
    cout << endl;

    cout << "old -> new" << endl;
    for (int old_id: old_ids) {
        int new_id = findNewId(old_list, new_list, old_id);
        cout << "'" << old_id << "': '" << new_id << "', ";
    }
}
