// g++ -O2 -march=native -I. bench_all.cpp -L. -lintegerb -o bench_all.exe
#include "integerb.h"
#include <chrono>
#include <cstdio>
#include <fstream>
#include <string>
#include <functional>

using clk = std::chrono::high_resolution_clock;
static double ms(clk::time_point a, clk::time_point b) {
    return std::chrono::duration<double, std::milli>(b - a).count();
}

int main() {
    floatb::fblen = 1000;

    auto mk = [](const char* s) {
        floatb x(s);
        return x;
    };

    struct Case {
        const char* name;
        floatb x;
        std::function<floatb(const floatb&)> f;
    };

    std::vector<Case> cases = {
        {"sin(0.2)",  mk("0.2"),  [](const floatb& x){ return fsin(x); }},
        {"cos(2.7)",  mk("2.7"),  [](const floatb& x){ return fcos(x); }},
        {"tan(5.6)",  mk("5.6"),  [](const floatb& x){ return ftan(x); }},
        {"exp(4.8)",  mk("4.8"),  [](const floatb& x){ return fexp(x); }},
        {"ln(123)",   mk("123"),  [](const floatb& x){ return fln(x); }},
    };

    std::ofstream out("result.txt");
    if (!out) { printf("cannot open result.txt\n"); return 1; }

    printf("%-12s %14s %10s\n", "case", "time(ms)", "validlen");
    for (auto& c : cases) {
        // 预热一次
        { volatile auto t = c.f(c.x); (void)t; }

        // 计时一次
        auto a = clk::now();
        floatb y = c.f(c.x);
        double t = ms(a, clk::now());

        printf("%-12s %14.3f %10d\n", c.name, t, y.validlength());

        // 写 result.txt：标题行 + 全精度数值行
        out << c.name << "\n" << y.tostring(0) << "\n";
    }
    out.close();
    printf("written result.txt\n");
    return 0;
}
/*


PS D:\mycpp> ./a
case               time(ms)   validlen
sin(0.2)             26.599       1000
cos(2.7)             25.831        999
tan(5.6)             26.200        999
exp(4.8)             17.078       1000
ln(123)              12.059       1000
written result.txt
PS D:\mycpp> & C:\Users\lenovo\AppData\Local\Programs\Python\Python39\python.exe d:/mycpp/test.py
=== 正确位数: good_limbs = -log(rel)/log(2^28), good_digits = -log10(rel) ===
case             good_limbs    good_digits
sin(0.2)            998.398         8415.3
cos(2.7)            997.578         8408.4
tan(5.6)            997.573         8408.4
exp(4.8)            997.817         8410.4
ln(123)             997.617         8408.8

=== MPFR 计时 (prec=28256 bits) ===
case               time(ms)
sin(0.2)             14.339
cos(2.7)             17.355
tan(5.6)             18.868
exp(4.8)             12.099
ln(123)               4.625



# & C:\Users\lenovo\AppData\Local\Programs\Python\Python39\python.exe d:/mycpp/test.py
# python check_lowprec.py —— 保留计时 + 正数 good_limbs/good_digits
import gmpy2
from gmpy2 import mpfr
import time

def good_limbs(a_cpp_str, a_std):
    """返回正数：和 MPFR 一致的正确 limb 数（越大越好）"""
    a = mpfr(a_cpp_str)
    diff = abs(a - a_std)
    if diff == 0:
        return float("inf")
    rel = diff / abs(a_std) if a_std != 0 else diff
    return float(-gmpy2.log(rel) / gmpy2.log(mpfr(2) ** 28))

def _f(name):
    if hasattr(gmpy2, name):
        return getattr(gmpy2, name)
    return lambda x: getattr(x, name)()

sin, cos, tan, exp, log = _f("sin"), _f("cos"), _f("tan"), _f("exp"), _f("log")

gmpy2.get_context().precision = 1000 * 28 + 256
gmpy2.get_context().round = gmpy2.RoundToNearest

std = {
    "sin(0.2)": sin(mpfr("0.2")),
    "cos(2.7)": cos(mpfr("2.7")),
    "tan(5.6)": tan(mpfr("5.6")),
    "exp(4.8)": exp(mpfr("4.8")),
    "ln(123)":  log(mpfr("123")),
}

with open(r"D:\mycpp\result.txt") as f:
    lines = [l.strip() for l in f if l.strip()]

# ---------- 1. 精度验证 ----------
print("=== 正确位数: good_limbs = -log(rel)/log(2^28), good_digits = -log10(rel) ===")
print("%-12s %14s %14s" % ("case", "good_limbs", "good_digits"))
for i in range(0, len(lines) - 1, 2):
    key = lines[i]
    if key not in std:
        print(f"{key}: 无标准值"); continue
    a_std = std[key]
    good = good_limbs(lines[i + 1], a_std)
    a = mpfr(lines[i + 1])
    diff = abs(a - a_std)
    rel = diff / abs(a_std) if a_std != 0 else diff
    good_digits = float(-gmpy2.log10(rel)) if rel != 0 else float("inf")
    print("%-12s %14.3f %14.1f" % (key, good, good_digits))

# ---------- 2. MPFR 计时（预热一次 + 计时一次） ----------
timing = [
    ("sin(0.2)",  sin, mpfr("0.2")),
    ("cos(2.7)",  cos, mpfr("2.7")),
    ("tan(5.6)",  tan, mpfr("5.6")),
    ("exp(4.8)",  exp, mpfr("4.8")),
    ("ln(123)",   log, mpfr("123")),
]

print("\n=== MPFR 计时 (prec=%d bits) ===" % (1000 * 28 + 256))
print("%-12s %14s" % ("case", "time(ms)"))
for name, fn, x in timing:
    _ = fn(x)                       # 预热一次
    t0 = time.perf_counter()
    y = fn(x)
    t = (time.perf_counter() - t0) * 1000.0
    print("%-12s %14.3f" % (name, t))

*/
