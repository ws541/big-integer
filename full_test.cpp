
#include "integer.h"
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
    integer r;
    integer c = a.divide(b, r);
    t2 = high_resolution_clock::now();
    ms = duration_cast<microseconds>(t2 - t1).count() / 1000.0;
    cout << "3. c = a / b = ";
    c.print();
    cout << "   (" << format_time(ms) << ")\n\n";
    // 4. 计算 d = sqroot(a)
    t1 = high_resolution_clock::now();
    integer d = sqroot(a);
    t2 = high_resolution_clock::now();
    ms = duration_cast<microseconds>(t2 - t1).count() / 1000.0;
    cout << "4. d = sqroot(a) = ";
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
