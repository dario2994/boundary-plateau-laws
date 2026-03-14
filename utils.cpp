#ifndef UTILS
#define UTILS

#define _USE_MATH_DEFINES
#include <bits/stdc++.h>
using namespace std;
typedef unsigned long long ULL;

//wrapper that checks if findIsomorphism returns a non-empty result.
template<typename T>
bool areIsomorphic(const T& A, const T& B) {
    return !findIsomorphism(A, B).empty();
}

//Standard matrix multiplication
vector<vector<ULL>> multiplyMatrices(const vector<vector<ULL>>& A,
                                     const vector<vector<ULL>>& B) {
    int n = A.size();
    assert(ssize(A[0]) == n and ssize(B) == n and ssize(B[0]) == n);
    vector<vector<ULL>> C(n, vector<ULL>(n, 0));
    for (int i = 0; i < n; i++) for (int j = 0; j < n; j++)
        for (int k = 0; k < n; k++) C[i][j] += A[i][k] * B[k][j];
    return C;
}

//Compute A^{256} for a matrix A
vector<vector<ULL>> largeMatrixPower(vector<vector<ULL>> mat) {
    for (int _ = 0; _ < 8; _++) mat = multiplyMatrices(mat, mat);
    return mat;
}

// Class to maintain a collection of objects without duplicates, when
// two isomorphic objects are considered equal.
// The underlying class T must have: isomorphismSignature and
// areIsomorphic.
template <typename T>
struct CollectionUpToIsomorphism {
    map<ULL, vector<T>> signature2list;

    void add(const T& x) {
        bool is_new = true;
        vector<T>& vec = signature2list[x.isomorphismSignature()];
        for (const T& y: vec) {
            if (areIsomorphic(x, y)) {
                is_new = false;
                break;
            }
        }
        if (is_new)
            vec.push_back(x);
        // The following assert is a safety measure to check that
        // the signature used is efficient.
        assert(ssize(vec) < 10);
    }

    vector<T> getList() {
        vector<T> res;
        for (const auto& pp: signature2list) for (const T& x: pp.second) {
            res.push_back(x);
        }
        return res;
    }
};


// Returns the time elapsed in nanoseconds from 1 January 1970, at 00:00:00.
long long get_time() {
    return chrono::duration_cast<chrono::nanoseconds>(
        chrono::steady_clock::now().time_since_epoch())
        .count();
}

// Tool for profiling code performance
class StopWatch {
  private:
    struct TimerData {
        long long total_time = 0;
        long long call_count = 0;
        bool is_active = false;
        bool is_anonymous = false;
    };
    inline static map<string,TimerData> timers;
    string name;

    static string randomName() {
        const string charset = "abcdefghijklmnopqrstuvwxyz";
        string result = "";
        for (int i = 0; i < 6; i++)
            result += charset[rand() % charset.size()];
        return result;
    }
  public:
    StopWatch(string _name): name(_name) {
        TimerData& td = timers[name];
        td.total_time -= get_time();
        td.call_count++;
        td.is_active = true;
    }

    StopWatch(): StopWatch(randomName()) {
        timers[name].is_anonymous = true;
    }

    void stop() {
        TimerData& td = timers[name];
        if (!td.is_active)
            return;
        td.is_active = false;
        td.total_time += get_time();
    }

    ~StopWatch() {
        stop();
    }

    static double get(string _name) { // Returns seconds
        assert(timers.count(_name));
        TimerData& td = timers[_name];
        long long elapsed = td.total_time;
        if (td.is_active) elapsed += get_time();
        double sec = elapsed / 1e9;
        return sec;
    }

    double get() {
        return get(name);
    }


    static void printState() {
        for (auto [name, td]: timers) {
            if (td.is_anonymous)
                continue;
            cout << name << ": " << get(name) << "s, "
                 << td.call_count << " calls" << endl;
        }
    }
};

//Error-handling tools
struct ErrorMessage {
    string error_str;
    ErrorMessage(string s): error_str(s) {}
    ErrorMessage(const stringstream& ss): error_str(ss.str()) {}
};

template <typename T>
class ErrorOr {
  private:
    bool is_error;
    string error_msg;
    unique_ptr<T> content;

    ErrorOr(bool err, std::string msg)
        : is_error(err), error_msg(std::move(msg)), content(nullptr) {}

  public:
    ErrorOr(T x)
        : is_error(false), error_msg(""), content(std::make_unique<T>(std::move(x))) {}

    ErrorOr(ErrorMessage EM): is_error(true), error_msg(std::move(EM.error_str)), content(nullptr) {}

    bool isOk() const { return !is_error; }

    T& getContent() const {
        assert(!is_error);
        return *content;
    }

    const string& getError() const {
        assert(is_error);
        return error_msg;
    }
};

template <typename T>
ostream& operator << (ostream& out, const ErrorOr<T>& error_or) {
    if (error_or.isOk())
        out << error_or.getContent();
    else
        out << error_or.getError();
    return out;
}

template <typename T>
T valueOrAssert(const ErrorOr<T>& error_or) {
    assert(error_or.isOk());
    return error_or.getContent();
}


#endif
