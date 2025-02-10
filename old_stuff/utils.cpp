#ifndef boundary_taylor_utils
#define boundary_taylor_utils

#include <bits/stdc++.h>
using namespace std;
typedef long long LL;

// Returns the time elapsed in nanoseconds from 1 January 1970, at 00:00:00.
LL get_time() {
    return chrono::duration_cast<chrono::nanoseconds>(
        chrono::steady_clock::now().time_since_epoch())
        .count();
}

// Returns a random integer from the range [l, r].
LL randint(LL l, LL r) {
	static mt19937_64 gen; // Deterministic on the local machine.
	return uniform_int_distribution<LL>(l, r)(gen);
}

// Debugging pretty print utilities.

template <typename Iterator>
string print_iterable(Iterator begin, Iterator end, string bef, string aft,
        function<string(typename iterator_traits<Iterator>::value_type)>* f = nullptr);

template <typename S, typename T>
ostream& operator <<(ostream& out, const pair<S, T>& p);

#define _define_print_container(Container, bef, aft)        \
template <typename T>                                       \
ostream& operator <<(ostream& out, const Container<T>& v) { \
    out << print_iterable(v.begin(), v.end(), bef, aft);    \
    return out;                                             \
}                                                           \

_define_print_container(vector, "[", "]");
_define_print_container(set, "{", "}");
_define_print_container(unordered_set, "{", "}");
_define_print_container(multiset, "{", "}");
_define_print_container(unordered_multiset, "{", "}");
_define_print_container(deque, "[", "]");

#define _define_print_dictionary(Dictionary, bef, aft)                      \
template <typename T1, typename T2>                                         \
ostream& operator <<(ostream& out, const Dictionary<T1, T2>& v) {           \
    out << print_iterable(v.begin(), v.end(), bef, aft);                    \
    return out;                                                             \
}                                                                           \

_define_print_dictionary(map, "{", "}");
_define_print_dictionary(unordered_map, "{", "}");

#define _define_print_container_adaptor(Adaptor, OP)                        \
template <typename T>                                                       \
ostream& operator <<(ostream& out, Adaptor<T> a) {                          \
    vector<T> v;                                                            \
    while (!a.empty()) v.push_back(a.OP()), a.pop();                        \
    out << v;                                                               \
    return out;                                                             \
}                                                                           \

_define_print_container_adaptor(stack, top)
_define_print_container_adaptor(queue, front)
_define_print_container_adaptor(priority_queue, top)

template <typename Iterator>
string print_iterable(Iterator begin, Iterator end, string bef, string aft,
        function<string(typename iterator_traits<Iterator>::value_type)>* f) {
    stringstream res;
    res << bef;
    bool is_first = true;
    for (Iterator it = begin; it != end; ++it) {
        if (!is_first) res << ", ";
        if (f != nullptr) res << (*f)(*it);
        else res << *it;
        is_first = false;
    }
    res << aft;
    return res.str();
}

template <typename S, typename T>
ostream& operator <<(ostream& out, const pair<S, T>& p) {
    out << "{" << p.first << ", " << p.second << "}";
    return out;
}

template<class TH>
void _dbg(const char* name, TH val){
    clog << name << ": " << val << endl;
}
template<class TH, class... TA>
void _dbg(const char* names, TH curr_val, TA... vals) {
    while(*names != ',') clog << *names++;
    clog << ": " << curr_val << ", ";
    _dbg(names+1, vals...);
}

#define dbg(...) _dbg(#__VA_ARGS__, __VA_ARGS__)


// Returns a vector u with the same length of v such that:
// 1. u[i] >= 0 for all i.
// 2. u[i] == u[j] <-> v[i] == v[j].
// 3. u is the lexicographically minimal vector with properties 1 and 2.
void minimal_permutation(vector<int>& v) {
    unordered_map<int,int> val2id;
    auto get_id = [&](int x) {
        if (val2id.count(x) == 0) {
            int sz = val2id.size();
            val2id[x] = sz;
        }
        return val2id[x];
    };
    for (int& x: v) x = get_id(x);
}

// Returns v[id] interpreting v in a circular fashion.
// id can be any integer.
template<typename T>
T& get(vector<T>& v, int id) {
    int sz = v.size();
    id %= sz;
    if (id < 0) id += sz;
    return v[id];
}

int get(const vector<int>& v, int id) {
    int sz = v.size();
    id %= sz;
    if (id < 0) id += sz;
    return v[id];
}

// This function removes duplicates, but instead of equality an isomorphism
// relation is considered.
// T is a type admitting a method T::isomorphism_class() returning a vector<int>
// which represents the isomorphism class of the instance of T.
template<typename T>
void remove_isomorphic_copies(vector<T>& vec) {
    set<vector<int>> isomorphism_classes;
    vector<T> res;
    for (T x: vec) {
        vector<int> ic = x.isomorphism_class();
        if (isomorphism_classes.count(ic)) continue;
        isomorphism_classes.insert(ic);
        res.push_back(x);
    }
    vec = res;
}

#endif
