#define _USE_MATH_DEFINES
#include <bits/stdc++.h>
using namespace std;
typedef unsigned long long ULL;

#include "derived_once_nets.cpp"


int findNewId(const vector<DerivedOnceNet>& old_list,
              const vector<DerivedOnceNet>& new_list,
              int old_id) {
    assert(0 <= old_id and old_id < ssize(old_list));
    const DerivedOnceNet& A = old_list[old_id];
    auto B = symmetrize(A);
    for (int new_id = 0; new_id < ssize(new_list); new_id++) {
        if (areIsomorphic(A, new_list[new_id])
            or areIsomorphic(B, new_list[new_id])) {
            return new_id;
        }
    }
    assert(0);
    return -1;
}

int main() {
    const vector<DerivedOnceNet>& old_list = parseDerivedOnceNets("data/old_derived_once_computer.txt");
    const vector<DerivedOnceNet>& new_list = parseDerivedOnceNets("data/derived_once_computer.txt");

    vector<int> old_ids = {
        72,     278,    434,    521,    691,
        936,    949,    1011,   1034,   1223,
        1351,   1405,   1514,   1829,   2017,
        2119,   2192,   2312,   2454,   2676,
        2700,   2769,   2799,   2891,   3048,
        3222,   3233,   3295,   3338};

    // vector<int> old_ids = {
    //     279, 19, 105, 1025, 115, 3091};

    vector<int> new_ids;

    cout << "old -> new" << endl;
    for (int old_id: old_ids) {
        int new_id = findNewId(old_list, new_list, old_id);
        cout << old_id << " -> " << new_id << endl;
        new_ids.push_back(new_id);
    }

    sort(new_ids.begin(), new_ids.end());
    cout << endl << "New ids list = " << endl;
    for (int new_id: new_ids)
        cout << new_id << ", ";
}
