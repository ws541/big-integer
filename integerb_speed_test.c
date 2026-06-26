#include "integerb.h"
#include <chrono>
#include <iostream>
#include <gmpxx.h>
#include <iomanip>

int main() {
    // ========== 1. Fibonacci Euclid 测试 ==========
    std::cout << "=== integer Euclid vs GMP exGCD (Fibonacci) ===" << std::endl;
    std::cout << "n\t\tsize\t\tinteger(us)\tGMP(us)\t\tGMP快几倍" << std::endl;
    std::cout << "-----------------------------------------------------------------" << std::endl;
    
    std::vector<int> fib_ns = {100, 500, 2500, 12500, 62500, 312500};
    std::vector<int> fib_int_times, fib_gmp_times;
    std::vector<int> fib_sizes;
    
    // integer Euclid
    for (int n : fib_ns) {
        pairs p = fib(n);
        fib_sizes.push_back(p.f.num.size());
        
        const int iterations = 4;
        integer s, t;
        auto start = std::chrono::high_resolution_clock::now();
        for (int i = 0; i < iterations; i++) {
            integer g = euclid(p.f, p.s, s, t);
        }
        auto end = std::chrono::high_resolution_clock::now();
        int int_time = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
        fib_int_times.push_back(int_time);
        if((p.f*s+p.s*t).num!=integer(1).num){std::cout<<"euclid bad";exit(0);}
    }
    
    // GMP exGCD
    for (int n : fib_ns) {
        mpz_class f, s;
        mpz_fib2_ui(f.get_mpz_t(), s.get_mpz_t(), n);
        
        const int iterations = 4;
        auto start = std::chrono::high_resolution_clock::now();
        for (int i = 0; i < iterations; i++) {
            mpz_class g, x, y;
            mpz_gcdext(g.get_mpz_t(), x.get_mpz_t(), y.get_mpz_t(), 
                       f.get_mpz_t(), s.get_mpz_t());
        }
        auto end = std::chrono::high_resolution_clock::now();
        int gmp_time = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
        fib_gmp_times.push_back(gmp_time);
    }
    
    for (size_t i = 0; i < fib_ns.size(); i++) {
        double ratio = (double)fib_int_times[i] / fib_gmp_times[i];
        std::cout << "n=" << fib_ns[i] 
                  << "\t\t" << fib_sizes[i]
                  << "\t\t" << fib_int_times[i]
                  << "\t\t" << fib_gmp_times[i]
                  << "\t\t" << std::fixed << std::setprecision(2) << ratio << "x" << std::endl;
    }
    
    // ========== 2. root7 测试 ==========
    std::cout << "\n=== integer root7 vs GMP root7 (fac n) ===" << std::endl;
    std::cout << "n\t\tsize\t\tinteger(us)\tGMP(us)\t\tGMP快几倍" << std::endl;
    std::cout << "-----------------------------------------------------------------" << std::endl;
    
    std::vector<int> root_ns = {100, 500, 2500, 12500, 62500};
    std::vector<int> root_int_times, root_gmp_times;
    std::vector<int> root_sizes;
    
    // integer root7
    for (int n : root_ns) {
        integer a = fac(n),r;
        root_sizes.push_back(a.num.size());
        
        const int iterations = 4;
        auto start = std::chrono::high_resolution_clock::now();
        for (int i = 0; i < iterations; i++) {
            r = a.root(7);
        }
        auto end = std::chrono::high_resolution_clock::now();
        int int_time = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
        root_int_times.push_back(int_time);
        integer c=power(r,7)-a,d=power(r+1,7)-a;
        if(c.sign>0&&c.num.back()||(d.sign<0||!d.num.back())){std::cout<<"root bad";exit(0);}
    }
    
    // GMP root7
    for (int n : root_ns) {
        mpz_class g;
        mpz_fac_ui(g.get_mpz_t(), n);
        
        const int iterations = 4;
        auto start = std::chrono::high_resolution_clock::now();
        for (int i = 0; i < iterations; i++) {
            mpz_class r;
            mpz_root(r.get_mpz_t(), g.get_mpz_t(), 7);
        }
        auto end = std::chrono::high_resolution_clock::now();
        int gmp_time = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
        root_gmp_times.push_back(gmp_time);
    }
    
    for (size_t i = 0; i < root_ns.size(); i++) {
        double ratio = (double)root_int_times[i] / root_gmp_times[i];
        std::cout << "n=" << root_ns[i] 
                  << "\t\t" << root_sizes[i]
                  << "\t\t" << root_int_times[i]
                  << "\t\t" << root_gmp_times[i]
                  << "\t\t" << std::fixed << std::setprecision(2) << ratio << "x" << std::endl;
    }
    
    // ========== 3. fac -> tostring -> construct 测试 ==========
    std::cout << "\n=== fac -> tostring -> construct (n=100000) ===" << std::endl;
    std::cout << "操作\t\tinteger(ms)\tGMP(ms)\t\tGMP快几倍" << std::endl;
    std::cout << "--------------------------------------------------------" << std::endl;
    
    // integer: fac(100000)
    auto start = std::chrono::high_resolution_clock::now();
    integer a = fac(100000);
    auto end = std::chrono::high_resolution_clock::now();
    long long int_fac = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    
    // GMP: fac(1000000)
    start = std::chrono::high_resolution_clock::now();
    mpz_class g;
    mpz_fac_ui(g.get_mpz_t(), 100000);
    end = std::chrono::high_resolution_clock::now();
    long long gmp_fac = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    double fac_ratio = (double)int_fac / gmp_fac;
    std::cout << "fac" << "\t\t" << int_fac << "\t\t" << gmp_fac << "\t\t" << std::fixed << std::setprecision(2) << fac_ratio << "x" << std::endl;
    
    // integer: tostring
    start = std::chrono::high_resolution_clock::now();
    std::string s = a.tostring();
    end = std::chrono::high_resolution_clock::now();
    long long int_str = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    
    // GMP: get_str
    start = std::chrono::high_resolution_clock::now();
    std::string gs = g.get_str();
    end = std::chrono::high_resolution_clock::now();
    long long gmp_str = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    double str_ratio = (double)int_str / gmp_str;
    std::cout << "tostring" << "\t" << int_str << "\t\t" << gmp_str << "\t\t" << std::fixed << std::setprecision(2) << str_ratio << "x" << std::endl;
    if(s!=gs){std::cout<<"tostring bad";exit(0);}
    // integer: construct from string
    start = std::chrono::high_resolution_clock::now();
    integer b(s);
    end = std::chrono::high_resolution_clock::now();
    long long int_construct = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    if(a.num!=b.num){std::cout<<"construct bad";exit(0);}
    // GMP: construct from string
    start = std::chrono::high_resolution_clock::now();
    mpz_class h(gs);
    end = std::chrono::high_resolution_clock::now();
    long long gmp_construct = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    double construct_ratio = (double)int_construct / gmp_construct;
    std::cout << "construct" << "\t" << int_construct << "\t\t" << gmp_construct << "\t\t" << std::fixed << std::setprecision(2) << construct_ratio << "x" << std::endl;
    
    // total
    long long int_total = int_fac + int_str + int_construct;
    long long gmp_total = gmp_fac + gmp_str + gmp_construct;
    double total_ratio = (double)int_total / gmp_total;
    std::cout << "--------------------------------------------------------" << std::endl;
    std::cout << "total" << "\t\t" << int_total << "\t\t" << gmp_total << "\t\t" << std::fixed << std::setprecision(2) << total_ratio << "x" << std::endl;
    
    return 0;
}
/*
root时间容易波动
PS D:\cpp\.vscode> ./a
=== integer Euclid vs GMP exGCD (Fibonacci) ===
n               size            integer(us)     GMP(us)         GMP快几倍
-----------------------------------------------------------------
n=100           3               45              65              0.69x
n=500           13              76              23              3.30x
n=2500          62              394             81              4.86x
n=12500         310             4262            803             5.31x
n=62500         1550            69480           14340           4.85x
n=312500                7749            491737          131442          3.74x

=== integer root7 vs GMP root7 (fac n) ===
n               size            integer(us)     GMP(us)         GMP快几倍
-----------------------------------------------------------------
n=100           19              228             23              9.91x
n=500           135             294             29              10.14x
n=2500          880             1724            142             12.14x
n=12500         5432            11226           1818            6.17x
n=62500         32342           74280           24993           2.97x

=== fac -> tostring -> construct (n=100000) ===
操作            integer(ms)     GMP(ms)         GMP快几倍
--------------------------------------------------------
fac             39              18              2.17x
tostring        219             67              3.27x
construct       59              30              1.97x
--------------------------------------------------------
total           317             115             2.76x
*/
