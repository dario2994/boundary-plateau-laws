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

    vector<int> interesting_ids = {};

    for (int old_id: interesting_ids) {
        int new_id = findNewId(old_list, new_list, old_id);
        cout << old_id << " -> " << new_id << endl;
    }
}
