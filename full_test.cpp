
#include "integerb.h"
#include <chrono>
#include <iomanip>
#include <sstream>
using namespace std;
using namespace chrono;
string format_time(double ms) {
    int sec = (int)(ms / 1000);
    int msec = (int)ms % 1000;
    stringstream ss;
    ss << sec << ".";
    if (msec < 10) ss << "00";
    else if (msec < 100) ss << "0";
    ss << msec << " s";
    return ss.str();
}
int main() {
    cout << "=====================================\n";
    cout << "   integer.h 测试结果\n";
    cout << "=====================================\n\n";
    auto total_start = high_resolution_clock::now();
    double ms;
    // 1. 计算 a = 1000000!
    auto t1 = high_resolution_clock::now();
    integer a = fac(1000000);
    auto t2 = high_resolution_clock::now();
    ms = duration_cast<microseconds>(t2 - t1).count() / 1000.0;
    cout << "1. a = 1000000! = ";
    a.print();
    cout << "   (" << format_time(ms) << ")\n\n";
    // 2. 计算 b = C(5000000, 2500000)
    t1 = high_resolution_clock::now();
    integer b = C(5000000, 2500000);
    t2 = high_resolution_clock::now();
    ms = duration_cast<microseconds>(t2 - t1).count() / 1000.0;
    cout << "2. b =C(5000000, 2500000) = ";
    b.print();
    cout << "   (" << format_time(ms) << ")\n\n";
    // 3. 计算 c = a / b
    t1 = high_resolution_clock::now();
    integer c = a/b;
    t2 = high_resolution_clock::now();
    ms = duration_cast<microseconds>(t2 - t1).count() / 1000.0;
    cout << "3. c = a / b = ";
    c.print();
    cout << "   (" << format_time(ms) << ")\n\n";
    // 4. 计算 d = sqrt(a)
    t1 = high_resolution_clock::now();
    integer d = a.fsqrt();
    t2 = high_resolution_clock::now();
    ms = duration_cast<microseconds>(t2 - t1).count() / 1000.0;
    cout << "4. d = sqrt(a) = ";
    d.print();
    cout << "   (" << format_time(ms) << ")\n\n";
    // 5. 计算 e = a 的 67 次方根
    t1 = high_resolution_clock::now();
    integer e67 = a.root(67);
    t2 = high_resolution_clock::now();
    ms = duration_cast<microseconds>(t2 - t1).count() / 1000.0;
    cout << "5. e = a 的 67 次方根 = ";
    e67.print();
    cout << "   (" << format_time(ms) << ")\n\n";
    // 6. 计算 f = gcd(100000!, 10^N)
    t1 = high_resolution_clock::now();
    integer f = fac(100000);
    integer g_shift = integer(1).shift(f.num.size());
    integer g = gcd(f, g_shift);
    t2 = high_resolution_clock::now();
    ms = duration_cast<microseconds>(t2 - t1).count() / 1000.0;
    cout << "6. g = gcd(100000!, 10^N) = ";
    g.print();
    cout << "   (" << format_time(ms) << ")\n\n";
    // 7. 计算 h = power(300!, 400!, 500!)
    t1 = high_resolution_clock::now();
    integer base = fac(300);
    integer exp = fac(400);
    integer mod = fac(500);
    integer h = power(base, exp, mod);
    t2 = high_resolution_clock::now();
    ms = duration_cast<microseconds>(t2 - t1).count() / 1000.0;
    cout << "7. h = power(300!, 400!, 500!) = ";
    h.print();
    cout << "   (" << format_time(ms) << ")\n\n";
    // 8. 找出 n! + 1 是质数的 n
    t1 = high_resolution_clock::now();
    cout << "8. n! + 1 是质数的 n (1~300)miller+lucas: ";
    integer fact = 1;
    for (int n = 1; n <= 300; n++) {
        fact = fact * integer(n);
        if (isprime(fact + 1)) {
            cout << n << " ";
        }
    }
    t2 = high_resolution_clock::now();
    ms = duration_cast<microseconds>(t2 - t1).count() / 1000.0;
    cout << "\n   (" << format_time(ms) << ")\n\n";
    // 9. 费马数 F7 的质因子
    t1 = high_resolution_clock::now();
    cout << "9. euler phi(F7 = 2^128 + 1)=";
    euler(power(2, 128) + 1).print();
    t2 = high_resolution_clock::now();
    ms = duration_cast<microseconds>(t2 - t1).count() / 1000.0;
    cout << "   (" << format_time(ms) << ")\n\n";
    auto total_end = high_resolution_clock::now();
    ms = duration_cast<microseconds>(total_end - total_start).count() / 1000.0;
    cout << "=====================================\n";
    cout << "总时间: " << format_time(ms) << "\n";
    cout << "=====================================\n";
    return 0;
}




#include"integer.h"
#include<chrono>
#include<iostream>
int main()
{
    integer n;
    n.get();//11111111111111111111111111122332231111111111111117791
    
    std::vector<integer> l;
    
    // Start timing
    auto start = std::chrono::high_resolution_clock::now();
    
    euler(n, &l);
    
    // End timing
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    
    // Print results
    for(auto& e : l) {
        e.print();
    }
    
    // Print timing information
    std::cout << "\nEuler function computation time: " 
              << duration.count() << " microseconds" << std::endl;
    std::cout << "Or: " << duration.count() / 1000.0 << " milliseconds" << std::endl;
    
    return 0;
}



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
    std::cout << "\n=== integer root3 vs GMP root3 (fac n) ===" << std::endl;
    std::cout << "n\t\tsize\t\tinteger(us)\tGMP(us)\t\tGMP快几倍" << std::endl;
    std::cout << "-----------------------------------------------------------------" << std::endl;
    
    std::vector<int> root_ns = {100, 500,1000,1500,2500, 12500, 62500};
    std::vector<int> root_int_times, root_gmp_times;
    std::vector<int> root_sizes;
    
    // integer root7
    for (int n : root_ns) {
        integer a = fac(n),r;
        root_sizes.push_back(a.num.size());
        
        const int iterations = 4;
        auto start = std::chrono::high_resolution_clock::now();
        for (int i = 0; i < iterations; i++) {
            r = a.root(3);
        }
        auto end = std::chrono::high_resolution_clock::now();
        int int_time = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
        root_int_times.push_back(int_time);
        integer c=power(r,3)-a,d=power(r+1,3)-a;
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
            mpz_root(r.get_mpz_t(), g.get_mpz_t(), 3);
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


#include "integerb.h"
#include <chrono>
#include <iostream>
int main()
{
    std::cout << "\n========== Index Calculus 离散对数求解 ==========\n" << std::endl;
    
    // 生成大质数 p (20! 附近的最小质数)
    integer p = fac(20);
    while(!isprime(p)) p.addsmall(1);
    
    // 底数 b = 13!
    integer b = fac(13);
    
    std::cout << "素数 p = " << p.tostring() << std::endl;
    std::cout << "底数 b = " << b.tostring() << std::endl;
    std::cout << "求解: log_g(b) mod p" << std::endl << std::endl;
    
    // ===== 阶段1: Index Calculus 预计算 =====
    std::cout << "[1] Index Calculus 预计算..." << std::endl;
    auto t1 = std::chrono::high_resolution_clock::now();
    integer g = primeroot(p);
    indexcalculus solver(p,g);
    auto t2 = std::chrono::high_resolution_clock::now();
    auto init_ms = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count() / 1000.0;
    std::cout << "    耗时: " << init_ms << " ms" << std::endl << std::endl;
    
    // ===== 阶段2: 求解离散对数 =====
    std::cout << "[2] 求解 log_g(b)..." << std::endl;
    auto t3 = std::chrono::high_resolution_clock::now();
    integer logb = solver.dlogr(b);
    auto t4 = std::chrono::high_resolution_clock::now();
    auto solve_ms = std::chrono::duration_cast<std::chrono::microseconds>(t4 - t3).count() / 1000.0;
    std::cout << "    log_g(b) = " << logb.tostring() << std::endl;
    std::cout << "    耗时: " << solve_ms << " ms" << std::endl << std::endl;
    
    // ===== 阶段3: 验证 =====
    std::cout << "[3] 验证: g^log ≡ b (mod p)..." << std::endl;
    auto t5 = std::chrono::high_resolution_clock::now();
    integer verify = power(g, logb, p);
    auto t6 = std::chrono::high_resolution_clock::now();
    auto verify_ms = std::chrono::duration_cast<std::chrono::microseconds>(t6 - t5).count() / 1000.0;
    
    std::cout << "    原根 g = " << g.tostring() << std::endl;
    std::cout << "    验证: " << (verify.num == b.num ? "✓ 通过" : "✗ 失败") << std::endl;
    std::cout << "    耗时: " << verify_ms << " ms" << std::endl << std::endl;
    
    // ===== 总结 =====
    std::cout << "========== 执行时间汇总 ==========" << std::endl;
    std::cout << "  Index Calculus 预计算: " << init_ms << " ms" << std::endl;
    std::cout << "  离散对数求解:          " << solve_ms << " ms" << std::endl;
    std::cout << "  验证:                  " << verify_ms << " ms" << std::endl;
    std::cout << "  总时间:                " << init_ms + solve_ms + verify_ms << " ms" << std::endl;
    std::cout << "===================================" << std::endl;
    
    std::cout << "\n结果: log_" << g.tostring() << "(" << b.tostring() << ") = " << logb.tostring() << std::endl;
    std::cout << "验证: " << g.tostring() << "^" << logb.tostring() << " ≡ " << b.tostring() << " (mod " << p.tostring() << ")" << std::endl;
    
    return 0;
}

