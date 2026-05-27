

/*
integer.h：支持快速乘法（karatsuba和FFT）、高效除法（burnikel-ziegler、牛顿迭代法和原生除法）、
组合数学函数（阶乘、排列、组合计算）、模运算系统（Montgomery模运算和通用模类）、高精度开平方、最大公约数（hgcd、lehmer包含扩展欧几里得算法）、
模逆元、快速模幂运算、素数判定（Miller-Rabin和Lucas测试）、jacobi、shanks模开方、质因数分解（Pollard Rho和二次筛法）、
欧拉函数、莫比乌斯函数计算以及原根求解,斐波那契数列。
只能作为演示使用，务必-O2或者release
*/
//[Console]::OutputEncoding = [System.Text.Encoding]::UTF8
//g++ a.cpp -O2 -I D:\mingw64download\mingw64\include -L D:\mingw64download\mingw64\lib -lgmp -lgmpxx -o a.exe
#ifndef _INTEGER_H_//-O2运行
#define _INTEGER_H_
#include<iostream>
#include<vector>
#include<cmath>
#include<complex>
typedef long long int ll;
typedef std::complex<double> cd;
constexpr int log2lenmax = 20, lenmax = 1 << log2lenmax, halflen = lenmax >> 1, Base = 100000000, Blen = 8, e4 = 10000;
//jacobi要求Blen>=3,quotient要求Blen>=2,sqroot要求Blen偶数,karamul不溢出要求Blen<=8,硬编码Blen=8
struct c2
{
private:
    struct d2//知乎学的rrii布局,ri压位卷积,循环分裂基写法
    {
    public:
        double a, b;
        d2() = default;
        constexpr d2(double x, double y) :a(x), b(y) {};
        d2 operator+(const d2& t) { return d2(a + t.a, b + t.b); }
        d2 operator-(const d2& t) { return d2(a - t.a, b - t.b); }
        d2 operator*(const d2& t) { return d2(a * t.a, b * t.b); }
        void  t() { double u = a + b, v = a - b; a = u, b = v; }
        d2 operator*(double t) { return d2(a * t, b * t); }
    };
    c2 operator*(const c2& t) { return c2(r * t.r, i * t.i); }
public:
    d2 r, i;
    c2() = default;
    constexpr c2(const d2& x, const d2& y) :r(x), i(y) {};
    c2(const cd& x, const cd& y) { r.a = real(x), i.a = imag(x), r.b = real(y), i.b = imag(y); }
    c2 operator+(const c2& t) { return c2(r + t.r, i + t.i); }
    c2 operator-(const c2& t) { return c2(r - t.r, i - t.i); }
    c2 mul(double t) { return c2(r * t, i * t); }
    c2 mul(const c2& t)
    {
        return c2(r * t.r - i * t.i, i * t.r + r * t.i);
    }
    c2 mulconj(const c2& t)
    {
        return c2(r * t.r + i * t.i, i * t.r - r * t.i);
    }
    void swapabconj()
    {
        std::swap(r.a, r.b);
        double u = -i.a;
        i.a = -i.b, i.b = u;
    }
    void bmulj() { double u = r.b; r.b = -i.b, i.b = u; }
    void bmulnegj() { double u = r.b; r.b = i.b, i.b = -u; }
    void mulj()
    {
        double u = r.a, v = r.b;
        r.a = -i.a, r.b = -i.b, i.a = u, i.b = v;
    }
    void mulnegj()
    {
        double u = r.a, v = r.b;
        r.a = i.a, r.b = i.b, i.a = -u, i.b = -v;
    }
    void t() { r.t(), i.t(); }
    void first(const c2& y)
    {
        double u = r.a * y.r.a + i.a * y.i.a;
        double v = r.a * y.i.a + i.a * y.r.a;
        r.a = u, i.a = v;
        u = r.b * y.r.b - i.b * y.i.b;
        v = r.b * y.i.b + i.b * y.r.b;
        r.b = u, i.b = v;
    }
    void second(const c2& y)
    {
        constexpr cd w = cd(0.7071067811865475, -0.7071067811865475), w2 = cd(0.3535533905932738, -0.3535533905932738);
        cd u = cd(r.a, i.a), v = cd(r.b, -i.b);
        cd xa = (u + v) * 0.5, xb = (u - v) * w2; xb = cd(std::imag(xb), -std::real(xb));
        u = cd(y.r.a, y.i.a), v = cd(y.r.b, -y.i.b);
        cd ya = (u + v) * 0.5, yb = (u - v) * w2; yb = cd(std::imag(yb), -std::real(yb));
        u = xa * ya + xb * yb, v = (xa * yb + xb * ya) * conj(w); v = cd(-std::imag(v), std::real(v));
        r.a = std::real(u) + std::real(v);
        i.a = std::imag(u) + std::imag(v);
        r.b = std::real(u) - std::real(v);
        i.b = std::imag(v) - std::imag(u);
    }
};
std::vector<double> table(int f)
{
    std::vector<double> l; l.resize(lenmax * 2);
    l[0] = l[2] = 1;
    for (int len = 4; len <= lenmax; len <<= 1)
    {
        int l2 = len >> 1;
        double t = (double)3.1415926535897932384626433832795 * f / l2, c0 = cos(t), s0 = sin(t);
        for (int i = l2, j = len; i < len; i += 2, j += 4)
        {
            double c = l[i], s = l[i + 1];
            l[j] = c;
            l[j + 1] = s;
            l[j + 2] = c * c0 - s * s0;
            l[j + 3] = c * s0 + c0 * s;
        }
    }
    for (int i = 0; i < l.size(); i += 4)
    {
        std::swap(l[i + 1], l[i + 2]);
    }
    return l;
}
std::vector<double>table1 = table(1), table3 = table(3);
void fft(c2* x, int n)
{
    for (int len = n; len > 2; len >>= 1)
    {
        int fftlen = len << 1, l = len >> 2;
        c2* t1 = reinterpret_cast<c2*>(&table1[fftlen]);
        c2* t3 = reinterpret_cast<c2*>(&table3[fftlen]);
        for (int i = len; i <= n; i <<= 2)
        {
            int gap = i << 1;
            for (int p = i - len; p < n; p += gap)
            {
                for (int k = 0; k < l; k++)
                {
                    int t = k + p;
                    c2 a = x[t], b = x[t + l], c = x[t + 2 * l], d = x[t + 3 * l];
                    c2 u = a - c, v = c2(b.i - d.i, d.r - b.r);
                    x[t] = a + c, x[t + l] = b + d, x[t + 2 * l] = (u + v).mulconj(t1[k]), x[t + 3 * l] = (u - v).mulconj(t3[k]);
                }
            }
        }
    }
    for (int i = 2; i <= n; i <<= 2)
    {
        int gap = i << 1;
        for (int t = i - 2; t < n; t += gap)
        {
            c2 u = x[t], v = x[t + 1];
            x[t] = u + v, x[t + 1] = u - v;
            x[t + 1].bmulnegj();
        }
    }
    for (int i = 0; i < n; i++) { x[i].t(); }
}
void ifft(c2* x, int n)
{
    for (int i = 0; i < n; i++) { x[i].t(); }
    for (int i = 2; i <= n; i <<= 2)
    {
        int gap = i << 1;
        for (int t = i - 2; t < n; t += gap)
        {
            c2 u = x[t], v = x[t + 1]; v.bmulj();
            x[t] = u + v, x[t + 1] = u - v;
        }
    }
    for (int len = 4; len <= n; len <<= 1)
    {
        int fftlen = len << 1, l = len >> 2;
        c2* t1 = reinterpret_cast<c2*>(&table1[fftlen]);
        c2* t3 = reinterpret_cast<c2*>(&table3[fftlen]);
        for (int i = len; i <= n; i <<= 2)
        {
            int gap = i << 1;
            for (int p = i - len; p < n; p += gap)
            {
                for (int k = 0; k < l; k++)
                {
                    int t = k + p;
                    c2 c = x[t + 2 * l].mul(t1[k]), d = x[t + 3 * l].mul(t3[k]);
                    c2 u = c + d, v = c2(d.i - c.i, c.r - d.r);
                    x[t + 2 * l] = x[t] - u, x[t + 3 * l] = x[t + l] - v;
                    x[t] = x[t] + u, x[t + l] = x[t + l] + v;
                }
            }
        }
    }
    double inv = 0.5L / n;
    for (int i = 0; i < n; i++) { x[i] = x[i].mul(inv); }
}
class integer
{
public:
    struct view
    {
    public:
        const int* ptr;
        int len;
        int sign;
        view(const integer& x) { ptr = x.num.data(), len = x.num.size(), sign = x.sign; }
        view(const int* a, int start, int end, int s) {
            while (end > start && a[end] == 0) { end--; }
            ptr = a + start;
            len = end - start + 1;
            sign = s;
        }
        view(const view& x, int leftshift) {//保证不没过长度
            ptr = x.ptr + leftshift;
            len = x.len - leftshift;
            if (len < 1 || leftshift < 0) { std::cout << "view"; exit(0); }
            sign = x.sign;
        }
    };
private:
    ll absll()const
    {
        ll a = 0, b = 1;
        for (ll e : num)
        {
            a += e * b;
            b *= Base;
        }
        return a;
    }
    static integer multiply(const int* a, int la, const int* b, int lb, int sign);
    int bz_init(const integer& b, integer& x, integer& y)const
    {
        int n = b.num.back(), f;
        if (n >= Base / 2)
        {
            f = 1;
            x.num = num, y.num = b.num;
        }
        else
        {
            f = Base / (n + 1);
            x = *this * f, y = b * f;
        }
        x.sign = y.sign = 1;
        return f;
    }
    void shiftadd(const view& b, int i)//不pushback1
    {
        int lb = b.len, k = 0;
        for (int j = 0; j < lb; j++)
        {
            num[j + i] += b.ptr[j] + k;
            if (num[i + j] >= Base) { num[i + j] -= Base; k = 1; }
            else { k = 0; }
        }
        if (k) { for (int j = i + lb; j < num.size() && ++num[j] == Base; j++) { num[j] = 0; } }
    }
    static  integer div_3n_2n(const view& a, const view& b, const view& b1, integer& r)
    {
        if (a.len < b.len + 30 || b.len < 30) { return div_native(a, b, r); }
        int n = b.len / 2, m = b1.len;
        integer q = div_4n_2n(view(a.ptr, n, a.len - 1, a.sign), b1, view(b1.ptr, m / 2, m - 1, b1.sign), r);
        int cnt = 0;
        integer product = karamul(q, b);
        r = a;
        int s = abssub(r.num.data(), r.num.size(), product.num.data(), product.num.size());
        while (s < 0)
        {
            if (cnt++ > 2)
            {
                std::cout << "div_3n_2n"; exit(0);
            }
            s = -abssub(r.num.data(), r.num.size(), b.ptr, b.len);
            q.addsmall(-1);
        }
        return q;
    }
    static  integer div_4n_2n(const view& a, const view& b, const view& b1, integer& r)
    {
        if (a.len < b.len + 30 || b.len < 30) return div_native(a, b, r);
        int n = b.len / 2;
        integer r32, q1 = div_3n_2n(view(a.ptr, n, a.len - 1, a.sign), b, b1, r32);
        r32.num.insert(r32.num.begin(), a.ptr, a.ptr + n);
        integer q; q.num.resize(a.len - b.len + 1);
        q.shiftadd(div_3n_2n(r32, b, b1, r), 0);
        if (q1.num.back()) {q.shiftadd(q1, n);}
        while (q.num.size() > 1 && !q.num.back()) q.num.pop_back();
        return q;
    }
    void quotient(const view& a, const view& b, integer& now)
    {
        int m = Base;
        int t = b.ptr[b.len - 1];
        while (t /= 10) { m /= 10; }
        if (!(m /= 10)) { m = 1; }
        ll m100 = m * 100, bb = m100 * b.ptr[b.len - 1] + (ll)m * b.ptr[b.len - 2] / (Base / 100);
        now.num.reserve(a.len + 1);
        now = a; now.num.push_back(0);
        for (int ln; (ln = now.num.size()) > 2 && (t = ln - b.len - 1) > -1; now.num.pop_back())
        {
            int q = (((ll)now.num.back() * Base + now.num[ln - 2]) * m100 + (ll)now.num[ln - 3] * m / (Base / 100)) / bb;
            if (!q) { continue; }
            integer prod = multiply(b.ptr, b.len, &q, 1, 1);
            int s = abssub(now.num.data() + t, b.len + 1, prod.num.data(), prod.num.size());
            int cnt = 0;
            while (s < 0)
            {
                s = -abssub(now.num.data() + t, b.len + 1, b.ptr, b.len);
                if (cnt++ > 2) { std::cout << "quotient"; exit(0); }
            }
            num[t] = q - cnt;
        }
        if (b.len == 2)
        {
            ll u = now.absll(), v = ((ll)b.ptr[1] * Base + b.ptr[0]);
            num[0] += u / v, now = u % v;
        }
        while (now.num.back() == 0 && now.num.size() > 1) { now.num.pop_back(); }
    }
    static integer div_native(const view& a, const view& b, integer& r)
    {
        if (!b.ptr[b.len - 1]) { std::cout << "div0"; exit(0); }
        integer c;
        c.num.resize(std::max(a.len - b.len + 1, 1));
        if (b.len == 1)//now*base不溢出
        {
            ll n = b.ptr[0], now = 0;
            for (int i = a.len - 1; i > -1; i--)
            {
                now = now * Base + a.ptr[i];
                c.num[i] = now / n;//(now*base+numi)/n<((n-1)*base+numi)/n<n*base/n=base
                now %= n;
            }
            r = now;
        }
        else
        {
            c.quotient(a, b, r);
        }
        while (c.num.size() > 1 && c.num.back() == 0) { c.num.pop_back(); }
        r.sign = a.sign;
        c.sign = a.sign * b.sign;
        return c;
    }
    static int abssub(int* a, int la, const int* b, int lb);
public:
    static integer addorsub(const int* a, int la, int asign, const int* b, int lb, int bsign, bool add);
    int sign = 1;
    std::vector<int> num;
    integer() {};
    integer(integer&& other)noexcept :sign(other.sign), num(std::move(other.num)) {}
    integer(const integer& other) :sign(other.sign), num(other.num) {}
    integer& operator=(integer&& other) noexcept
    {
        sign = other.sign;
        num = std::move(other.num);
        return *this;
    }
    integer& operator=(const integer& other)
    {
        sign = other.sign;
        num = other.num;
        return *this;
    }
    integer(const view& v)
    {
        num.assign(v.ptr, v.ptr + v.len);
        sign = v.sign;
    }
    integer(ll n)
    {
        if (n < 0) { n = -n; sign = -1; }
        do
        {
            num.push_back(n % Base);
            n /= Base;
        } while (n);
    }
    integer(const std::string& s)
    {
        int i = s.size() - Blen;
        for (; i > 0; i -= Blen)
        {
            int a = s[i] - 48;
            for (int j = 1; j < Blen; j++)
            {
                a *= 10;
                a += s[i + j] - 48;
            }
            num.push_back(a);
        }
        i += Blen;
        int a = 0, j = 0;
        if (s[0] == '-') { j = 1; sign = -1; }
        for (; j < i; j++)
        {
            a *= 10;
            a += s[j] - 48;
        }
        if (a || num.empty()) { num.push_back(a); }
    }
    ll toll()const
    {
        if (absbigger(integer(9223372036854775807), 0)) { std::cout << "toll"; exit(0); }
        return sign * absll();
    }
    void get()
    {
        std::string s;
        std::cin >> s;
        *this = integer(s);
    }
    void print(int n = 10)const
    {
        if (sign != 1) { std::cout << "-"; }
        int l = num.size(), cnt = 0; int a[Blen];
        a[Blen - 1] = num.back();
        for (int j = Blen - 1; j > 0; j--)
        {
            a[j - 1] = a[j] / 10;
        }
        for (; cnt < Blen - 1 && a[cnt] == 0; cnt++);
        if (l > n) { std::cout << a[cnt++] << "."; }
        for (int j = cnt; j < Blen; j++)
        {
            std::cout << a[j] % 10;
        }
        for (int i = 1; i < n && i < l; i++)
        {
            a[Blen - 1] = num[l - 1 - i];
            for (int j = Blen - 1; j > 0; j--)
            {
                a[j - 1] = a[j] / 10;
            }
            for (int j = 0; j < Blen; j++)
            {
                std::cout << a[j] % 10;
            }
        }
        if (l > n)
        {
            std::cout << "e" << l * Blen - cnt;
        }
        std::cout << std::endl;
    }
    bool absbigger(const integer& that, bool equal)const
    {
        int n = num.size();
        if (n == that.num.size())
        {
            int i = n - 1;
            for (; i > -1 && num[i] == that.num[i]; i--);
            if (i == -1) { return equal; }
            return num[i] > that.num[i];
        }
        return n > that.num.size();
    }
    integer operator-()const
    {
        integer a = *this; a.sign = -sign;
        return a;
    }
    integer operator+(const integer& that)const
    {
        return addorsub(num.data(), num.size(), sign, that.num.data(), that.num.size(), that.sign, 1);
    }
    integer operator-(const integer& that)const
    {
        return addorsub(num.data(), num.size(), sign, that.num.data(), that.num.size(), that.sign, 0);
    }
    static integer karamul(view  a, view b);
    static integer fftmul(const view& a, const view& b);
    static double* e4rrii(const view& a, int n);
    integer operator*(const integer& that)const
    {
        return karamul(*this, that);
    }
    void reciprocal(integer& xt, int l) const  // xt=base^t/this,,t最终=l,对商长度和*this长度有最小值要求
    {
        int k = 1;
        int t = num.size() + k;
        int a = num.size(), b;
        integer tmp;
        xt = div_native(integer(1).shift(3), view(num.data(), a - 2, a - 1, sign), tmp);
        int l1 = 0;
        if (l > 200)
        {
            l1 = l - num.size();
            while (l1 > 40) { l1 = (l1 + 1) / 2; }l1 += num.size();
        }
        while (t < l)
        {
            if (t + k >= l) { k = l - t; }  // 最终微调
            else if (l1 && k + t >= l1) { k = l1 - t; l1 = 0; }//中间级微调,使得最终xt几乎充分利用
            tmp = karamul(xt, xt);//想要触发fft平方优化
            int need = xt.num.size() + k + 2;  // 向(xt * 2).shift(k)约xt.size+k对齐
            b = std::max(a - need, 0);
            int b1 = std::max((int)tmp.num.size() - need, 0);
            tmp = karamul(view(*this, b), view(tmp, b1));
            xt = (xt + xt).shift(k);
            view tmpview(tmp, -(b + b1 + k - t));
            abssub(xt.num.data(), xt.num.size(), tmpview.ptr, tmpview.len);
            while (xt.num.back() == 0 && xt.num.size() > 1) { xt.num.pop_back(); }
            t += k;  // 新增k个待计算（有效数字增长一倍)
            k = xt.num.size() - 1;
        }
        xt.addsmall(-sign);  // 不能大于,认为最多差1
    }
    integer fsqrt() const
    {
        int k = 1;
        int ns = num.size();
        int t = ns + 2 * k;
        int l = 2 * (ns + 2);
        int b = std::max(ns - (k + 2), 0);
        integer xt = (ll)sqrt((ll)Base * Base / num.back()) * 10000, tmp = shift(-b), y;
        while ((y = (xt * 3 - ((tmp * xt) * (xt * xt)).shift(-t + b)) / 2).num != xt.num)
        {
            std::swap(xt, y);
        }
        int l1 = 0;
        if (l > 200)
        {
            l1 = l - ns;
            while (l1 > 40) { l1 = (l1 + 1) / 2; }l1 += ns;
        }
        while (t < l)
        {
            if (t + 2 * k > l) { k = (l - t + 1) / 2; }
            else if (l1 && t + 2 * k > l1) { k = (l1 - t + 1) / 2; l1 = 0; }
            tmp = karamul(xt, xt);
            int need = k + xt.num.size() + 2;
            b = std::max(ns - need, 0);
            int b1 = std::max((int)tmp.num.size() - need, 0);
            tmp = karamul(view(*this, b), view(tmp, b1));
            int b2 = std::max((int)tmp.num.size() - need, 0);
            tmp = karamul(xt, view(tmp, b2));
            view tmpview(tmp, -(k - t + b + b1 + b2));//为了保证是左移需要*this足够大,小的用除法版本
            xt = (xt * 3).shift(k);
            abssub(xt.num.data(), xt.num.size(), tmpview.ptr, tmpview.len);
            while (xt.num.back() == 0 && xt.num.size() > 1) { xt.num.pop_back(); }
            xt = xt / 2;
            t += 2 * k;
            k = xt.num.size() - 1;
        }
        b = num.size() - xt.num.size() - 2; if (b < 0) { b = 0; }
        integer r = karamul(xt, view(*this, b)).shift(b - t / 2);
        if (r.num.size() != ns / 2) { r = r / 10000; }
        integer y2 = r * r - *this;
        if (y2.num.size() > ns / 2 + 1) { std::cout << "fsqrt"; exit(0); }
        if (y2.sign == 1)
        {
            if (y2.num.back()) { r.addsmall(-1); }
            return r;
        }
        y2 = y2 + r * 2;
        if (y2.sign == -1 && y2.num.back()) { r.addsmall(1); }
        return r;
    }
    integer div_newton(const integer& that, integer& r)const//余数符号同被除数(由图像)
    {
        r.num.clear();
        int la = num.size(), lb = that.num.size();
        int n = std::min(24, la / lb);
        int l = (la - lb) / n + 3;
        while (l >= 512) { l >>= 1; }n += l < 285 && l>255;
        l = (la - lb) / n + lb + 1;//(l-that.size)*n=size-that.size,保护+1;
        integer xt;
        that.reciprocal(xt, l);
        integer res(*this), q;
        q.num.assign(num.size() - that.num.size() + 1, 0);
        for (int i = 0, b = that.num.size() - 2; i < n && res.num.size()>that.num.size() + 4; i++)
        {
            int lr = res.num.size(), start = lr - l;
            if (start < 0) { start = 0; }
            view cutview(res.num.data(), start, lr - 1, res.sign);
            integer p = karamul(xt, view(cutview, b));
            view pview(p, l - b);
            integer now = karamul(that, pview);
            if (abssub(res.num.data() + start, std::min(l, lr), now.num.data(), now.num.size()) < 0) { std::cout << "div_newton"; exit(0); }
            q.shiftadd(pview, start);
            while (!res.num.back() && res.num.size() > 1) { res.num.pop_back(); }
        }
        if (res.num.size() > that.num.size() + 4) { std::cout << "div_newton"; exit(0); }
        q.shiftadd(div_native(res, that, r), 0);
        while (q.num.size() > 1 && q.num.back() == 0) { q.num.pop_back(); }
        q.sign = sign * that.sign;
        return q;
    }
    integer divide(const integer& that, integer& r)const
    {
        int a = num.size(), b = that.num.size(), c = a - b;
        if (a > 50 && b > 1 && c > 15 && (c > 70 || b * c > 1600))
        {
            return div_newton(that, r);
        }
        return div_native(*this, that, r);
    }
    integer div_bz(const integer& b, integer& r) const//不如牛顿迭代法快,同样有最小长度要求
    {
        integer y;
        int v = bz_init(b, r, y);
        int m = y.num.size(), n = 2 * m;
        integer q, r1;
        view y1(y.num.data(), m / 2, m - 1, y.sign);
        q.num.assign(r.num.size() - y.num.size() + 1, 0); m = 1;
        while (m && r.num.size() > y.num.size())
        {
            int lr = r.num.size(); m = lr - n;
            if (m < 0) { m = 0; }
            view cutview(r.num.data(), m, lr - 1, 1);
            q.shiftadd(div_4n_2n(cutview, y, y1, r1), m);
            for (int j = 0; j < r1.num.size(); j++)
            {
                r.num[m + j] = r1.num[j];
            }
            r.num.resize(m + r1.num.size());
            while (!r.num.back() && r.num.size() > 1) { r.num.pop_back(); }
        }
        while (q.num.size() > 1 && q.num.back() == 0) { q.num.pop_back(); }
        q.sign = sign * b.sign;
        r = r / v;
        r.sign = sign;
        return q;
    }
    integer operator/(const integer& that)const
    {
        integer r;
        return divide(that, r);
    }
    integer operator%(const integer& that)const
    {
        integer r;
        divide(that, r);
        return r;
    }
    integer shift(int n)const
    {
        if (!num.back() || !n) { return *this; }
        if ((int)num.size() + n < 1) { return 0; }
        integer c;
        if (n > 0)
        {
            c.num.resize(n, 0);
            c.num.insert(c.num.end(), num.begin(), num.end());
        }
        else
        {
            c.num.insert(c.num.begin(), num.begin() - n, num.end());
        }
        c.sign = sign;
        return c;
    }
    integer cut(int start, int end)const
    {
        integer a;
        while (end > start && num[end] == 0) { end--; }
        a.num.assign(num.begin() + start, num.begin() + end + 1);
        a.sign = sign;
        return a;
    }
    void addsmall(int n)//n绝对值不超过Base
    {
        num[0] += sign == 1 ? n : -n;
        int l = num.size();
        if (num[0] >= Base)
        {
            num[0] -= Base;
            int i;
            for (i = 1; i < l && ++num[i] == Base; i++) { num[i] = 0; }
            if (i == l) { num.push_back(1); }
        }
        else if (num[0] < 0)
        {
            if (l == 1) { num[0] *= -1; sign *= -1; return; }
            num[0] += Base;
            int i;
            for (i = 1; i < l && --num[i] == -1; i++) { num[i] = Base - 1; }
            if (!num.back()) { num.pop_back(); }
        }
    }
    static void mul_core(std::vector<integer>& level, int& zeros)
    {
        int s, s2;
        while ((s = level.size()) > 1)
        {
            s2 = s / 2;
            std::vector<integer> nextlevel;
            for (int k = 0; k < s2; k++)
            {
                integer& a = level[k];
                integer& b = level[k + s2];
                int i = 0;
                while (i < (int)a.num.size() - 1 && a.num[i] == 0) i++;
                zeros += i; view a_view(a, i);
                i = 0;
                while (i < (int)b.num.size() - 1 && b.num[i] == 0) i++;
                zeros += i; view b_view(b, i);
                nextlevel.push_back(karamul(a_view, b_view));
            }
            if (s % 2) { nextlevel.push_back(level.back()); }
            level = std::move(nextlevel);
        }
    }
};
int integer::abssub(int* a, int la, const int* b, int lb)
{
    if (la < lb) { std::cout << "abssub"; exit(0); }
    for (int i = 0; i < lb; i++)
    {
        a[i] -= b[i];
    }
    while (a[la - 1] == 0 && la > 1) { la--; }
    int s; bool k = 0;
    if (a[la - 1] < 0)
    {
        s = -1;
        for (int i = 0; i < la; i++)
        {
            a[i] = -a[i] - k;
            k = a[i] < 0;
            a[i] += k * Base;
        }
    }
    else
    {
        int i = 0; s = 1;
        for (; i < lb && i < la; i++)
        {
            a[i] -= k;
            k = a[i] < 0;
            a[i] += k * Base;
        }
        for (; i < la && k; i++)
        {
            a[i] -= k;
            k = a[i] < 0;
            a[i] += k * Base;
        }
    }
    return s;
}
integer integer::addorsub(const int* a, int la, int asign, const int* b, int lb, int bsign, bool add)
{
    int change = 1;
    if (la < lb) { std::swap(a, b); std::swap(la, lb); std::swap(asign, bsign); change = -1; }
    bool k = 0; int j;
    integer result;
    const bool absadd = (asign == bsign) == add;
    result.num.reserve(la + absadd);
    result.num.assign(a, a + la);
    result.sign = asign;
    if (absadd)
    {
        for (int i = 0; i < lb; i++)
        {
            result.num[i] += b[i] + k;
            k = result.num[i] >= Base;
            result.num[i] -= k * Base;
        }
        const int lr = result.num.size();
        for (j = lb; j < lr && k; j++)
        {
            result.num[j] += k;
            k = result.num[j] >= Base;
            result.num[j] -= k * Base;
        }
        if (k) { result.num.push_back(1); }
    }
    else
    {
        result.sign *= abssub(result.num.data(), result.num.size(), b, lb);
        while (result.num.back() == 0 && result.num.size() > 1) { result.num.pop_back(); }
    }
    if (!add) { result.sign *= change; }
    return result;
}
integer integer::multiply(const int* a, int la, const int* b, int lb, int sign)//保证la不小于lb
{
    if ((la + a[0] < 2) || (lb + b[0] < 2)) { return 0; }
    int n = la - 1, m = lb - 1, l = n + m + 1;
    integer c;
    c.sign = sign;
    c.num.resize(l + 1);
    ll k = 0;
    int i;
    for (i = 0; i < m; i++)
    {
        for (int j = i; j >= 0; j--)
        {
            k += (ll)a[i - j] * b[j];
        }
        c.num[i] = k % Base;
        k = k / Base;
    }
    for (; i < n; i++)
    {
        for (int j = m; j >= 0; j--)
        {
            k += (ll)a[i - j] * b[j];
        }
        c.num[i] = k % Base;
        k = k / Base;
    }
    for (; i < l; i++)
    {
        for (int j = m; j >= i - n; j--)
        {
            k += (ll)a[i - j] * b[j];
        }
        c.num[i] = k % Base;
        k = k / Base;
    }
    c.num[l] = k;
    while (c.num.size() > 1 && c.num.back() == 0) { c.num.pop_back(); }
    return c;
}
integer integer::karamul(view a, view b)
{
    if (a.len < b.len) { std::swap(a, b); }
    if (b.len < 80) return multiply(a.ptr, a.len, b.ptr, b.len, a.sign * b.sign);
    if (a.len + b.len <= halflen && b.len > 160) return fftmul(a, b);
    int n = (a.len + 1) / 2;
    view a1(a.ptr, n, a.len - 1, 1);
    view a0(a.ptr, 0, n - 1, 1);
    integer ans;
    if (b.len <= n)
    {
        view bb = b; bb.sign = 1;
        ans = karamul(a0, bb);
        ans.num.resize(a.len + b.len);
        ans.shiftadd(karamul(a1, bb), n);
    }
    else
    {
        view b1(b.ptr, n, b.len - 1, 1);
        view b0(b.ptr, 0, n - 1, 1);
        integer c1 = karamul(a1, b1);
        integer c0 = karamul(a0, b0);

        integer a0_a1 = addorsub(a0.ptr, a0.len, 1, a1.ptr, a1.len, 1, 0);
        integer b0_b1 = addorsub(b0.ptr, b0.len, 1, b1.ptr, b1.len, 1, 0);

        ans.num = c0.num;
        ans.num.resize(a.len + b.len);
        ans.shiftadd(c0 + c1 - karamul(a0_a1, b0_b1), n);
        ans.shiftadd(c1, 2 * n);
    }
    while (ans.num.size() > 1 && ans.num.back() == 0) ans.num.pop_back();
    ans.sign = a.sign * b.sign;
    return ans;
}
double* integer::e4rrii(const view& a, int n)
{
    double* y = new double[n]();
    int i = 0;
    for (; i + 1 < a.len; i += 2)
    {
        y[2 * i] = a.ptr[i] % e4;
        y[2 * i + 2] = a.ptr[i] / e4;
        y[2 * i + 1] = a.ptr[i + 1] % e4;
        y[2 * i + 3] = a.ptr[i + 1] / e4;
    }
    if (i != a.len)
    {
        y[2 * i] = a.ptr[i] % e4;
        y[2 * i + 2] = a.ptr[i] / e4;
    }
    return y;
}
std::vector<cd> unit = []() {
    std::vector<cd> l;
    for (int i = 0; i < log2lenmax; i++) {
        double t = -(double)3.1415926535897932384626433832795 / (1 << (i + 1));
        l.push_back(cd(cos(t), sin(t)));
    }
    return l;
    }();
integer integer::fftmul(const  view& a, const view& b)
{
    int la = a.len, lb = b.len;
    int m = (la + lb) << 1, n = 256;
    while (n < m) { n <<= 1; }
    if (n > lenmax || n < 0) { std::cout << "fftmul"; exit(0); }
    int n2 = n >> 1, n4 = n >> 2;
    double* x0 = e4rrii(a, n), * y0;
    c2* xx = reinterpret_cast<c2*>(&x0[0]), * yy;
    fft(xx, n4);
    const bool same = a.ptr == b.ptr && a.len == b.len;
    if (same) { y0 = x0, yy = xx; }
    else { y0 = e4rrii(b, n), yy = reinterpret_cast<c2*>(&y0[0]); fft(yy, n4); }
    xx[0].first(yy[0]); xx[1].second(yy[1]);
    int l = 0; m = n2;
    while (!(m & 1)) { m >>= 1, l++; }
    int lk = l - 3;
    std::vector<cd> tb; tb.assign(log2lenmax, 0); int index;
    auto next = [&tb, &index, &n4](int& k)->cd {
        int bit = n4, cnt = 0;
        while (k & bit) { k ^= bit, bit >>= 1; cnt++; }k ^= bit;
        index -= cnt; tb[index + 1] = tb[index] * unit[cnt];
        return tb[++index];
        };
    for (int len = 2; len < n4; len <<= 1, lk--)
    {
        int k = 1 << lk; index = 0, tb[0] = unit[l - 1 - lk];
        c2 w(tb[0], next(k));
        for (int s = len, e = (len << 1) - 1;;)
        {
            c2 w2 = w.mul(0.5);
            c2 u = xx[s], v = xx[e]; v.swapabconj();
            c2 xa = (u + v).mul(0.5), xb = (u - v).mul(w2); xb.mulj();
            u = yy[s], v = yy[e]; v.swapabconj();
            c2 ya = (u + v).mul(0.5), yb = (u - v).mul(w2); yb.mulj();
            u = xa.mul(ya) + xb.mul(yb), v = (xa.mul(yb) + xb.mul(ya)).mulconj(w); v.mulnegj();
            xx[s] = u + v, xx[e] = u - v; xx[e].swapabconj();
            if (--e > ++s)
            {
                cd t = next(k); w = c2(t, next(k));
            }
            else { break; }
        }
    }if (!same) { delete[]y0; }
    ifft(xx, n4);
    integer c;
    c.num.reserve(m = la + lb); m = (m + 1) / 2;
    c.sign = a.sign * b.sign;
    ll k = 0;
    for (int i = 0; i < m; i++)
    {
        ll t = (ll)(xx[i].r.a + 0.5) + k;
        ll pre = t % e4;
        k = t / e4;
        t = (ll)(xx[i].i.a + 0.5) + k;
        k = t / e4;
        c.num.push_back((t % e4) * e4 + pre);
        t = (ll)(xx[i].r.b + 0.5) + k;
        pre = t % e4;
        k = t / e4;
        t = (ll)(xx[i].i.b + 0.5) + k;
        k = t / e4;
        c.num.push_back((t % e4) * e4 + pre);
    }
    delete[]x0;
    while (c.num.size() > 1 && c.num.back() == 0) c.num.pop_back();
    return c;
}
integer power(const integer& b, const integer& c)
{
    if (c.sign == -1 || c.num.size() > 1) { std::cout << "power"; exit(0); }
    int n = c.num[0];
    integer a(b), ans = n & 1 ? b : 1;
    while (n >>= 1)
    {
        a = a * a;
        if (n & 1) { ans = ans * a; }
    }
    return ans;
}
integer C_core(int e, int s, int  e_s, int& zeros)//s<=e
{
    int a, b = (e + e % 2) / 2;
    std::vector<int>u(b, 1);
    for (int i = 1; 2 * (a = 2 * i * (i + 1)) <= e - 1; i++)
    {
        if (u[i])
        {
            for (int j = a; j < b; j += 2 * i + 1) { u[j] = 0; }
        }
    }u[2] = 0;
    std::vector<int>p;
    for (int i = 1; i < b; i++) { if (u[i]) { p.push_back(2 * i + 1); } }
    auto cnt = [](int n, int p)->int {
        int a = 0; ll c = p;
        while (n >= c) { a += n / (int)c; c *= p; }
        return a;
        };
    a = cnt(e, p[0]) - cnt(s, p[0]) - cnt(e_s, p[0]);
    u.clear();
    for (int i = 1; i < p.size(); i++)
    {
        b = cnt(e, p[i]) - cnt(s, p[i]) - cnt(e_s, p[i]);
        if (a != b)
        {
            u.push_back(i); u.push_back(a); a = b;
        }
    }u.push_back(p.size()), u.push_back(a);
    auto mul_lambda = [&p](int s, int e)->integer {
        std::vector<integer>l;
        for (e--; s < e; s += 2)
        {
            l.push_back((ll)p[s] * p[s + 1]);
        }
        if (s == e) { l.push_back(p[e]); }
        integer::mul_core(l, s);
        return l.back();
        };
    std::vector<integer>l;
    a = cnt(e, 2) - cnt(s, 2) - cnt(e_s, 2), b = cnt(e, 5) - cnt(s, 5) - cnt(e_s, 5);
    if (a - b < 0) { l.push_back(power(5, b - a)); b = a; }
    else { l.push_back(power(2, a - b)); }
    for (int i = 0, j = 0; i < u.size() - 2; j = u[i], i += 2)
    {
        l.push_back(power(mul_lambda(j, u[i]), u[i + 1]));
    }
    integer::mul_core(l, s);
    integer r = (l.back() * power(10, b % Blen));
    if (u.size() > 2) { r = r * power(mul_lambda(u[u.size() - 4], u[u.size() - 2]), u[u.size() - 1]); }
    zeros += b / Blen;
    return r;
}
integer mul(std::vector<integer>& level)
{
    int zeros = 0; integer::mul_core(level, zeros); return level.back().shift(zeros);
}
integer A_core(int e, int s, int& zeros)//s<=e
{
    if (s > 10000 && (e >> 2) < s) { return C_core(e, e - s, 0, zeros); }
    std::vector<integer> level;
    for (s = e - s + 1; s < e; s += 2)
    {
        level.push_back((ll)s * (s + 1));
    }
    if (s == e) { level.push_back(e); }
    integer::mul_core(level, zeros);
    return level.back();
}
integer A(int e, int s)//s<=e
{
    if (s < 1) { return 1; }int zeros = 0;
    return A_core(e, s, zeros).shift(zeros);
}
integer fac(int n)
{
    if (n < 10000) { return A(n, n - 1); }
    int zeros = 0, n2 = n / 2;
    integer f = A_core(n2, n2 - 1, zeros); f = f * f; zeros <<= 1;
    if (n & 1) { f = f * n; n--; }
    return (C_core(n, n2, n - n2, zeros) * f).shift(zeros);
}
integer C(int e, int s)//s<=e
{
    s = std::min(s, e - s);
    if (s < 2) { return s ? e : 1; }
    if (s * 100 < e * 3) { return A(e, s) / A(s, s - 1); }
    int zeros = 0; return  C_core(e, s, e - s, zeros).shift(zeros);
}
void gcdshift(integer& x, integer& y, const integer& q)
{
    std::swap(x, y), y = y - q * x;
}
void gcdshift1(integer& x, integer& y, const integer& A, const integer& B, const integer& C, const integer& D)
{
    integer tmp;
    std::swap(tmp, x), x = A * tmp + C * y, y = B * tmp + D * y;
}
void gcdshift2(integer& a, integer& b, integer& c, integer& d,
    const integer& A, const integer& B, const integer& C, const integer& D)
{
    integer u, v;
    std::swap(u, a); a = u * A + b * C; b = u * B + b * D;
    std::swap(v, c); c = v * A + d * C; d = v * B + d * D;
}
bool lehmer(integer& x, integer& y, integer* a, integer* b, integer* c, integer* d)
{
    if (a) { a->sign *= x.sign; c->sign *= x.sign; b->sign *= y.sign; d->sign *= y.sign; }
    x.sign = y.sign = 1;
    ll xh = (ll)x.num.back() * Base + x.num[x.num.size() - 2], yh = y.num.back();
    if (x.num.size() == y.num.size() && y.num.size() > 1) { yh = yh * Base + y.num[y.num.size() - 2]; }
    ll A = 1, B = 0, C = 0, D = 1, q, u, v;
    while ((u = yh + B) && (v = yh + D) && (q = (xh + A) / u) == (xh + C) / v)
    {
        u = A, A = B, B = u - B * q;
        u = C, C = D, D = u - C * q;
        u = xh, xh = yh, yh = u - yh * q;
    }
    if (B)
    {
        gcdshift1(x, y, A, B, C, D);
        if (a) { gcdshift2(*a, *b, *c, *d, A, B, C, D); }
    }
    return B;
}
void hgcd(const integer::view& x, const integer::view& y, integer& a, integer& b, integer& c, integer& d)//多项式版本粗暴映射
{
    int m = x.len / 2;
    if (m < 50)
    {
        integer xn(x), yn(y);
        while (yn.num.back() && yn.num.size() >= m)
        {
            int l = yn.num.size();
            if (!(l > m && xn.num.size() < l + 2 && lehmer(xn, yn, &a, &b, &c, &d)))
            {
                integer r, q = xn.divide(yn, r);
                gcdshift(a, b, q);
                gcdshift(c, d, q);
                std::swap(xn, yn); std::swap(yn, r);
            }
        }return;
    }
    else if (y.len - 1 < m) { return; }
    integer A(1), B(0), C(0), D(1), xn(x), yn(y);
    hgcd(integer::view(xn.num.data(), m, xn.num.size() - 1, xn.sign), integer::view(yn.num.data(), m, yn.num.size() - 1, yn.sign), A, B, C, D);
    gcdshift1(xn, yn, A, B, C, D);
    gcdshift2(a, b, c, d, A, B, C, D);
    if (!yn.num.back() || yn.num.size() - 1 < m) { return; }
    integer r;
    A = 0, B = C = 1, D = -xn.divide(yn, r);
    if (yn.num.size() >= 2 * m) { std::cout << "hgcd"; exit(0); }
    if (r.num.back() && r.num.size() > m)
    {
        int k = 2 * m - yn.num.size();
        hgcd(integer::view(yn.num.data(), k, yn.num.size() - 1, yn.sign), integer::view(r.num.data(), k, r.num.size() - 1, r.sign), A, B, C, D);
    }
    gcdshift2(a, b, c, d, A, B, C, D);
}
integer gcd(const integer& m, const integer& n)
{
    integer x(m), y(n);
    if (y.absbigger(x, 0)) { std::swap(x, y); }
    while (y.num.back())
    {
        while (y.num.size() > 200 && 10 * y.num.size() > 9 * x.num.size())
        {
            integer A(1), B(0), C(0), D(1);
            hgcd(x, y, A, B, C, D);
            gcdshift1(x, y, A, B, C, D);
        }
        if (y.num.back())
        {
            int l = x.num.size();
            if (!(l > 1 && l < y.num.size() + 2 && lehmer(x, y, 0, 0, 0, 0)))
            {
                integer r, q1 = x.divide(y, r);
                std::swap(x, y); std::swap(y, r);
            }
        }
    }
    x.sign = 1;
    return x;
}
integer euclid(const integer& m, const integer& n, integer& a, integer& c)
{
    integer x(m), y(n);
    bool change = 0;
    if (y.absbigger(x, 0)) { std::swap(x, y); change = 1; }
    a = 1, c = 0;
    integer b(0), d(1);
    while (y.num.size() > 200)
    {
        while (10 * y.num.size() > 9 * x.num.size())
        {
            integer A(1), B(0), C(0), D(1);
            hgcd(x, y, A, B, C, D);
            gcdshift1(x, y, A, B, C, D);
            gcdshift2(a, b, c, d, A, B, C, D);
        }
        if (y.num.size() > 200)
        {
            integer r, q = x.divide(y, r);
            gcdshift(a, b, q);
            gcdshift(c, d, q);
            std::swap(x, y); std::swap(y, r);
        }
    }
    if (y.num.back())
    {
        integer A(1), B(0), C(0), D(1);
        while (y.num.back())
        {
            int l = x.num.size();
            if (!(l > 1 && l < y.num.size() + 2 && lehmer(x, y, &A, &B, &C, &D)))
            {
                integer r, q = x.divide(y, r);
                gcdshift(A, B, q);
                gcdshift(C, D, q);
                std::swap(x, y); std::swap(y, r);
            }
        }
        gcdshift2(a, b, c, d, A, B, C, D);
    }
    if (x.sign == -1)
    {
        x.sign = 1;
        a.sign = -a.sign;
        c.sign = -c.sign;
    }
    if (change) { std::swap(a, c); }
    return x;
}
integer inv(const integer& a, const integer& p)//p>0
{
    integer s, t, g = euclid(p, a, s, t);
    if (!(g.num.size() == 1 && g.num[0] == 1) || !t.num.back()) { return 0; }
    if (t.sign == 1) { return t % p; }
    return t % p + p;
}
integer sqroot(const integer& n, bool trust = 1)
{
    if (n.sign == -1 && n.num.back()) { std::cout << "sqroot"; exit(0); }
    if (n.num.size() < 3)
    {
        ll a = n.toll(), b = sqrt(a);
        if (n.num.size() > 1 && b * b > a) { b -= 1; }//(Base-1)^2-1
        return b;
    }
    if (trust && n.num.size() > 30) { return n.fsqrt(); }
    integer y, x;
    int a; ll t = sqrt((ll)n.num.back() * Base + n.num[n.num.size() - 2]) + 1;
    if (n.num.size() > 6)//trust=0除法版本牛顿迭代慢
    {
        a = n.num.size();
        int k = 1, b = a - 2 * k, l = (a + 2) / 2;
        x = t;
        while ((k = x.num.size()) < l)
        {
            int j = 2 * k > l ? l - k : k;
            x = x.shift(j);
            b = a - 2 * (k + j);//b会小于0所以不能view,因为必须从大于一侧逼近
            x = (x + n.shift(-b) / x) / 2;
        }
        b = (n.num.size() - 1) * Blen, a = n.num.back();
        while (a = a / 10) { b++; }b = b / 2;
        a = x.num.back(), k = (x.num.size() - 1) * Blen;
        while (a = a / 10) { k++; }
        x = x / power(10, k - b);
        a = 2;
    }
    else
    {
        if (n.num.size() & 1) { t *= 10000; }
        x = integer(t).shift(n.num.size() / 2 - 1);
        a = -1;
    }
    while (x.absbigger(y = (x + n / x) / 2, 0))
    {
        if (a-- == 0) { std::cout << "sqroot"; exit(0); }
        std::swap(x, y);
    }
    return x;
}
class mont
{
private:
    integer a, b;
    bool normal;
    int m;
    bool flag(const integer& x)
    {
        if (x.num.size() == m)
        {
            for (auto& e : x.num)
            {
                if (e < Base - 1) { return 0; }
            }
            return 1;
        }
        return x.num.size() > m;
    }
    void sub(integer& x)
    {
        int i = 0;
        for (; i < m && x.num[i] == Base - 1; i++)
        {
            x.num[i] = 0;
        }
        if (i != m)
        {
            x.num[i]++;
            for (i = m; i < x.num.size() - 1 && x.num[i] == 0; i++) { x.num[i] = Base - 1; }
            x.num[i] -= 1;
        }
        while (x.num.back() == 0 && x.num.size() > 1) { x.num.pop_back(); }
    }
    integer fastmod(const integer& x)
    {
        integer ans = x.cut(0, std::min((int)x.num.size() - 1, m - 1));
        if (normal) { return ans; }
        if (m < x.num.size()) {
            integer::view tmp(x.num.data(), m, x.num.size() - 1, x.sign);
            ans = integer::addorsub(ans.num.data(), ans.num.size(), ans.sign, tmp.ptr, tmp.len, tmp.sign, 1);
        }
        int cnt = 0;
        while (flag(ans))
        {
            if (cnt++ > 2) { std::cout << "fastmod"; exit(0); }
            sub(ans);
        }
        return ans;
    }
    integer fastdiv(const integer& x)
    {
        integer ans;
        if (m > x.num.size() - 1) { ans = 0; }
        else { ans = x.cut(m, x.num.size() - 1); }
        if (normal) { return ans; }
        integer remain;
        if (ans.num.back())
        {
            integer::view tmp(x.num.data(), 0, m - 1, x.sign);
            remain = integer::addorsub(ans.num.data(), ans.num.size(), ans.sign, tmp.ptr, tmp.len, tmp.sign, 1);
        }
        else { remain = x; }
        int cnt = 0;
        while (flag(remain))
        {
            if (cnt++ > 2) { std::cout << "fastmod"; exit(0); }
            sub(remain);
        }
        ans.addsmall(cnt);
        return ans;
    }
public:
    integer f, p;
    mont() {}
    void init(bool type)
    {
        normal = type;
        m = p.num.size();
        integer n;
        if (!normal) { n.num.assign(m, Base - 1); }
        else { n = integer(1).shift(m); }
        a = n - inv(p, n);
        f = n % p;
        b = f * f % p;
    }
    integer out(const integer& x)
    {
        integer t = fastdiv(x + p * fastmod(fastmod(x) * a));
        if (t.absbigger(p, 1)) { return t - p; }
        return t;
    }
    integer in(const integer& x)
    {
        integer tmp = x * b % p;
        return out(tmp.sign > 0 || tmp.num.back() == 0 ? tmp : tmp + p);
    }
    integer inv_inside(const integer& x)
    {
        return in(inv(out(x), p));
    }
    integer pow_binary(const integer& x, const std::vector<int>& l)
    {
        if (l.size() == 1 && l[0] == 0) { return f; }
        integer u(x), v;
        bool k = 1;
        const int ls = l.size() - 1;
        for (int i = 0; i <= ls; i++)
        {
            ll m = l[i];
            for (int j = 0; j < 26; j++)
            {
                if (m & 1)
                {
                    if (k) { v = u; k = 0; }
                    else { v = out(u * v); }
                }
                if (!(m >>= 1) && i == ls) { break; }
                u = out(u * u);
            }
        }
        return v;
    }
    integer add(const integer& u, const integer& v)
    {
        integer ans(u + v);
        if (ans.sign == -1 && ans.num.back()) { ans = ans + p; }
        else if (ans.absbigger(p, 1)) { ans = ans - p; }
        return ans;
    }
    integer sub(const integer& u, const integer& v)
    {
        integer ans(u - v);
        if (ans.sign == -1 && ans.num.back()) { ans = ans + p; }
        else if (ans.absbigger(p, 1)) { ans = ans - p; }
        return ans;
    }
    integer pow_exponent(const integer& u, const std::vector<int>& l, int w = 5)
    {
        if (l.empty()) { return f; }
        std::vector<integer>b; b.push_back(u);
        integer a = out(u * u); w = (1 << (w - 1)) - 1;
        for (int i = 0; i < w; i++)
        {
            b.push_back(out(a * b.back()));
        }
        const int ls = l.size();
        a = b[l[ls - 2]];
        for (int j = 0; j < l[ls - 1]; j++) { a = out(a * a); }
        for (int i = ls - 4; i > -1; i -= 2)
        {
            a = out(a * b[l[i]]);
            for (int j = 0; j < l[i + 1]; j++)
            {
                a = out(a * a);
            }
        }
        return a;
    }
};
struct pairs
{
    integer f;
    integer s;
};
class mod
{
private:
    integer key;
    bool a1, b1;
    mont b;
public:
    mont a;
    mod(const integer& x, int step = 128)//>1
    {
        int m = x.num[0];
        if (m % 2 && m % 5)
        {
            a.p = x;
            a.init(1); a1 = 1, b1 = 0;
            return;
        }
        int i = 0;
        for (; i < x.num.size() - 1 && x.num[i] == 0; i++);
        m = x.num[i];
        int j = 1;
        while (m % 10 == 0) { m /= 10; j *= 10; }
        m = m % 5 ? 2 : 5;
        integer p = power(m, step), res(j), q = x.cut(i, x.num.size() - 1) / res, g;
        do
        {
            p = p * p;
            g = gcd(p, q);
        } while (g.num == p.num);
        p = q / g;
        if (!(p.num.size() == 1 && p.num[0] == 1)) { a.p = std::move(p); a.init(1); a1 = 1; }
        else { a1 = 0; }
        p = (g * res).shift(i);
        if (!(p.num.size() == 1 && p.num[0] == 1)) { b.p = std::move(p); b.init(0); b1 = 1; }
        else { b1 = 0; }
        if (a1 && b1) { key = inv(a.p, b.p); }
    }
    pairs in(const integer& x)
    {
        pairs ans;
        if (a1)
        {
            ans.f = a.in(x % a.p);
        }
        if (b1)
        {
            ans.s = b.in(x % b.p);
        }
        return ans;
    }
    pairs mul(const pairs& u, const pairs& v)
    {
        pairs ans;
        if (a1)
        {
            ans.f = a.out(u.f * v.f);
        }
        if (b1)
        {
            ans.s = b.out(u.s * v.s);
        }
        return ans;
    }
    pairs add(const pairs& u, const pairs& v)
    {
        pairs ans;
        if (a1)
        {
            ans.f = a.add(u.f, v.f);
        }
        if (b1)
        {
            ans.s = b.add(u.s, v.s);
        }
        return ans;
    }
    pairs sub(const pairs& u, const pairs& v)
    {
        pairs ans;
        if (a1)
        {
            ans.f = a.sub(u.f, v.f);
        }
        if (b1)
        {
            ans.s = b.sub(u.s, v.s);
        }
        return ans;
    }
    integer out(const pairs& u)
    {
        if (!a1) { return  b.out(u.s); }
        integer x = a.out(u.f);
        if (!b1) { return x; }
        integer tmp = (b.out(u.s) - x) * key % b.p;
        if (tmp.num.back() && tmp.sign == -1) { tmp = tmp + b.p; }
        return x + a.p * tmp;
    }
    pairs one()
    {
        pairs ans;
        if (a1)
        {
            ans.f = a.f;
        }
        if (b1)
        {
            ans.s = b.f;
        }
        return ans;
    }
    pairs zero()
    {
        pairs ans;
        if (a1)
        {
            ans.f = 0;
        }
        if (b1)
        {
            ans.s = 0;
        }
        return ans;
    }
    pairs inv_inside(const pairs& x)
    {
        pairs ans;
        if (a1)
        {
            ans.f = a.inv_inside(x.f);
            if (!ans.f.num.back()) { return zero(); }
        }
        if (b1)
        {
            ans.s = b.inv_inside(x.s);
            if (!ans.s.num.back()) { return zero(); }
        }
        return ans;
    }
    bool iszero(const pairs& x)
    {
        return  (!(a1 && x.f.num.back())) && !(b1 && x.s.num.back());
    }
    pairs pow_binary(const pairs& x, const std::vector<int>& l)
    {
        pairs ans;
        if (a1)
        {
            ans.f = a.pow_binary(x.f, l);
        }
        if (b1)
        {
            ans.s = b.pow_binary(x.s, l);
        }
        return ans;
    }
    pairs pow_exponent(const pairs& x, const std::vector<int>& l, int w = 5)
    {
        pairs ans;
        if (a1)
        {
            ans.f = a.pow_exponent(x.f, l, w);
        }
        if (b1)
        {
            ans.s = b.pow_exponent(x.s, l, w);
        }
        return ans;
    }
};
std::vector<int>tobinary(const integer& n)
{
    integer b(1 << 26), r, a = n.divide(b, r);
    std::vector<int>l; l.push_back(r.num[0]);
    while (a.num.back())
    {
        a = a.divide(b, r);
        l.push_back(r.num[0]);
    }
    return l;
}
std::vector<int> exponent(integer a, int w = 5)
{
    integer c(1 << 26), r;
    int step = 0, now = 0, cnt = 0;
    std::vector<int>l;
    while (a.num.back())
    {
        a = a.divide(c, r);
        int m = r.num[0];
        for (int i = 0; i < 26 && (a.num.back() || m); i++, m >>= 1)
        {
            int x = m & 1;
            if (x || now)
            {
                now += x << cnt;
                if (++cnt == w)
                {
                    l.push_back(now / 2);
                    l.push_back(step);
                    cnt = now = 0; step = w;
                }
            }
            else { step++; }
        }
    }
    if (now)
    {
        l.push_back(now / 2);
        l.push_back(step);
    }
    return l;
}
integer power(const integer& x, const integer& n, const integer& p)//p>1
{
    int m = n.num.back();
    if (n.num.size() == 1)
    {
        integer a(x % p), ans = m & 1 ? a : 1;
        while (m >>= 1)
        {
            a = a * a % p;
            if (m & 1) { ans = ans * a % p; }
        }
        return ans.sign == -1 && ans.num.back() ? ans + p : ans;
    }
    mod q(p); int w = std::min((int)n.num.size(), 4) + 1;
    return q.out(q.pow_exponent(q.in(x), exponent(n, w), w));
}
int ctz(integer& d)//d改变
{
    integer x(d), b(1 << 26), r;
    int s = 0, j = 0;
    int m;
    while (x.num.back())
    {
        x = x.divide(b, r);
        if (m = r.num[0]) { break; }
        s++;
    }
    while (!(m & 1)) { m >>= 1; j++; }
    d = (x * b + r) / (1 << j);
    return s * 26 + j;
}
bool inside(int k, const std::initializer_list<int>& t) {
    for (int v : t) {
        if (v == k) { return 1; }
    }
    return 0;
}
bool miller(const integer& n, const std::initializer_list<int> t = { 2, 3, 5, 7, 11, 13, 17 })
{
    if (n.sign == -1 || n.num.size() == 1 && n.num[0] == 1) { return 0; }
    if (n.num.size() == 1 && inside(n.num[0], t)) { return 1; }
    if (n.num[0] % 2 == 0 || n.num[0] % 5 == 0) { return 0; }
    integer d(n - 1);
    int s = ctz(d);
    mont q; q.p = n; q.init(1);
    integer flag = q.p - q.f;
    int w = n.num.size();
    w = 1 + (w > 1) + (w > 2) + (w > 3) + (w > 4);
    std::vector<int> l = exponent(d, w);
    for (int e : t)
    {
        integer now = q.pow_exponent(q.in(e), l, w);
        if (now.num == q.f.num || now.num == flag.num) { continue; }
        int k;
        for (k = 0; k < s; k++)
        {
            now = q.out(now * now);
            if (now.num == flag.num) { break; }
        }
        if (k == s) { return 0; }
    }
    return 1;
}
int jacobi(const integer& a, const integer& b)
{
    if (b.num[0] % 2 == 0 || b.sign == -1) { std::cout << "jacobi"; exit(0); }
    integer x(a % b), y(b);
    if (x.sign == -1 && x.num.back()) { x = x + y; }
    int m = 1;
    while (x.num.back())
    {
        int j = y.num[0] % 8;
        if (x.num[0] % 2 == 0)
        {
            int k = ctz(x);
            int s = (j == 1 || j == 7) ? 1 : -1;
            if (s == -1 && k % 2) { m = -m; }
        }
        if (x.num.size() == 1 && x.num[0] == 1) { return m; }
        int i = x.num[0] % 4;
        if (i == 3 && j % 4 == 3) { m = -m; }
        std::swap(x, y);
        x = x % y;
    }
    return 0;
}
bool lucas(const integer& x)
{
    if (x.num.size() == 1 && (x.num[0] == 2 || x.num[0] == 5)) { return 1; }
    integer u;
    if (x.num[0] % 2 == 0 || x.num[0] % 5 == 0 || ((u = sqroot(x)) * u).num == x.num) { return 0; }
    int t = 5, s = 1;
    while (jacobi(t * s, x) != -1)
    {
        t += 2; s = -s;
    }
    std::vector<int>l = exponent(x + 1, 1);
    mont q; q.p = x; q.init(1); u = q.f;
    integer v = u, d = q.in(t * s), qk = q.in((1 - t * s) / 4), q0 = qk, inv2 = q.in((x + 1) / 2);
    const int ls = l.size();
    for (int j = 0; j < l[ls - 1]; j++)
    {
        u = q.out(u * v);
        v = q.sub(q.out(v * v), q.add(qk, qk));
        qk = q.out(qk * qk);
    }
    for (int i = ls - 4; i > -1; i -= 2)
    {
        integer uu; std::swap(uu, u);
        u = q.out(inv2 * q.add(uu, v));
        v = q.out(q.add(q.out(d * uu), v) * inv2);
        qk = q.out(qk * q0);
        for (int j = 0; j < l[i + 1]; j++)
        {
            u = q.out(u * v);
            v = q.sub(q.out(v * v), q.add(qk, qk));
            qk = q.out(qk * qk);
        }
    }
    return !u.num.back();
}
bool isprime(const integer& x)
{
    if (x.absbigger(341550071728321, 1)) { return miller(x, { 2,3 }) && lucas(x); }
    return miller(x);
}
integer shanks_sqrt(const integer& a, const integer& p)//p奇质数
{
    if (p.num.size() == 1 && p.num[0] == 5)
    {
        int m = (a.num[0] * a.sign) % 5;
        if (m < 0) { m += 5; }
        return (m == 1) + 2 * (m == 4);
    }
    int n = 2;
    while (jacobi(n, p) != -1) { n++; }
    integer t = p - 1;
    int e = ctz(t);
    mont q; q.p = p; q.init(1);
    integer y = q.pow_binary(q.in(n), tobinary(t)), a0 = q.in(a),
        x = q.pow_binary(a0, tobinary((t - 1) / 2)), b = q.out(a0 * q.out(x * x));
    x = q.out(a0 * x);
    while (b.num != q.f.num)
    {
        int m = 0;
        integer k = b;
        do
        {
            if (m == e - 1) { return 0; }
            k = q.out(k * k); m++;
        } while (k.num != q.f.num);
        k = y;
        for (int i = m + 1; i < e; i++)
        {
            k = q.out(k * k);
        }
        y = q.out(k * k), e = m, x = q.out(x * k), b = q.out(b * y);
    }
    return q.out(x);
}
integer factor(const integer& n, bool& pollard)//不质数幂
{
    if (n.num[0] % 2 == 0) { return 2; }
    if (n.num[0] % 5 == 0) { return 5; }
    integer g;
    auto check_lambda = [&g, &n]() -> bool {
        return (g.num.size() > 1 || g.num[0] > 1) && g.num != n.num;
        };
    int B = 1, ok = 0;
    if (pollard)
    {
        do
        {
            integer a(1), b(1), c(B), p(1);
            ll cnt = 0; bool go = 1;
            for (int step = 1; go; step <<= 1, a = b)
            {
                for (int i = 0; i < step && go; i++)
                {
                    b = (b * b + c) % n;
                    if (b.num == a.num)
                    {
                        g = gcd(p, n);
                        go = false; break;
                    }
                    p = p * (a - b) % n;
                    if (!p.num.back()) { return gcd(a - b, n); }
                    if (++cnt % 50 == 0)
                    {
                        g = gcd(p, n);
                        if (check_lambda()) { return g; }
                        else if (ok++ > 1500) { go = false; }
                    }
                }
            }
            B++;
        } while (ok < 1501 && !check_lambda());
        if (check_lambda()) { return g; }
        pollard = 0;
    }
    if (n.num.size() < 3) { std::cout << "pollard"; exit(0); }
    integer d = sqroot(n);
    mont q;
    std::vector<int>prime;
    std::vector<int>root;
    std::vector<float>logp;
    struct equa
    {
        integer x, y;
        std::vector<unsigned long long>mask;
    public:
        equa() {};
        equa(int ns, int ps) { x.num.resize(ns), y.num.resize(ns), mask.resize(ps / 64 + (bool)(ps % 64 != 0)), x.sign = -1; }
    };
    std::vector<equa>smooth;//qs知乎学的,e45之内保证不超过10s,超过e47时间指数增长
    std::vector<equa>flaw; int fs = 0;
    integer X, Y; int ta = 0, tb = 0, tc = 0; float bound;
    auto insert_lambda =
        [&q, &prime, &smooth, &flaw, &fs, &X, &Y, &n, &ok, &ta, &tb, &tc, &bound](const integer& x, integer y) -> integer {
        int  ps = prime.size();
        equa e; e.y = y; e.mask.assign(ps / 64 + (bool)(ps % 64 != 0), 0);
        for (int i = 0; i < ps && (y.num.size() > 1 || y.num[0] > 1); i++)
        {
            integer r, p = prime[i], q = y.divide(p, r);
            int cnt = 0;
            while (!r.num.back())
            {
                cnt++;
                std::swap(y, q);
                q = y.divide(p, r);
            }
            e.mask[i / 64] |= ((unsigned long long)(cnt & 1)) << (i % 64);
        }
        if (y.num.size() > 1 || y.num[0] > Base / 20) { tb++; return 0; }
        int j = 0;
        if (y.num[0] > 1)
        {
            tb++; for (; j < fs && flaw[j].x.sign != y.num[0]; j++); if (j >= flaw.size() || !miller(y.num[0], { 2 })) { return 0; }
        }
        else { ta++; }
        e.x = q.in(x);
        if (y.num[0] > 1) { e.y = e.y / y; }
        int  c = 0, r, cend = e.mask.size() * 64;
        for (; c < cend; c += 64)
        {
            unsigned long long ei = e.mask[c / 64];
            for (r = 0; r < ps - c && r < 64; r++, ei >>= 1)
            {
                if (ei & 1) { e.y = e.y / prime[c + r]; }
            }
        }
        e.y.sign = 1, e.y = sqroot(e.y) % n, e.y.sign = y.sign;
        if (y.num[0] > 1)
        {
            equa& s = flaw[j];
            if (s.x.sign < 0) { e.x.sign = y.num[0]; flaw[j] = e; fs++; return 0; }
            tc++;
            s.x.sign = 1; e.x = q.out(e.x * s.x); s.x.sign = y.num[0];
            e.y = e.y * (s.y * y.num[0]);
            for (c = 0; c < cend; c += 64)
            {
                unsigned long long ei = e.mask[c / 64], si = s.mask[c / 64], u = ei & si;
                for (r = 0; r < ps - c && r < 64; r++, u >>= 1)
                {
                    if (u & 1) { e.y = e.y * prime[c + r]; }
                }
                e.mask[c / 64] ^= si;
            }
            e.y = e.y % n;
        }
        for (c = 0; c < cend && !e.mask[c / 64]; c += 64);
        unsigned long long bit;
        if (c < cend) { for (bit = e.mask[c / 64], r = 0; r < ps - c && r < 64 && !(bit & 1); r++, bit >>= 1); }
        int index;
        while ((index = c + r) < ps && smooth[index].x.sign > 0)
        {
            equa& s = smooth[index];
            e.x = q.out(e.x * s.x);
            int pre = c;
            std::vector<integer>l;
            for (; c < cend; c += 64)
            {
                unsigned long long ei = e.mask[c / 64], si = s.mask[c / 64], u = ei & si;
                for (r = 0; r < ps - c && r < 64; r++, u >>= 1)
                {
                    if (u & 1) { l.push_back(prime[c + r]); }
                }
                e.mask[c / 64] ^= si;
            }
            e.y = mul(l) % n * (e.y * s.y) % n;
            for (c = pre; c < cend && !e.mask[c / 64]; c += 64);
            if (c < cend) { for (bit = e.mask[c / 64], r = 0; r < ps - c && r < 64 && !(bit & 1); r++, bit >>= 1); }
        }
        if (index >= ps)
        {
            if (e.y.sign < 0)
            {
                if (X.num.empty()) { X = e.x, Y = e.y; return 0; }
                e.x = q.out(e.x * X), e.y = e.y * Y;
            }ok--;
            return gcd(q.out(e.x) + e.y, n);
        }
        smooth[index] = e;
        //std::cout << "\r    relation:" << ta + tc << " need:" << prime.size() << " |from " <<tb << " fail extract " << fs << " flaw and " << tc << " relation\r" << std::flush;
        const int level = 4;
        if (tb > (ta << (level + 1))) { bound += 0.1; }if ((ta << level) > tb) { bound -= 0.1; }//倾向于错误多些,快
        return 0;
        };
    auto getlog = [](const integer& x)->float {
        int l = x.num.size();
        return l == 1 ? log(x.num[0]) : log((ll)x.num.back() * Base + x.num[l - 2]) + (l - 2) * log(Base);
        };
    B = pow(getlog(d) - 2, 1.7 + 0.1 * n.num.size());
    root.assign((B + B % 2) / 2, 1);
    for (int i = 1, a; 2 * (a = 2 * i * (i + 1)) <= B - 1; i++)
    {
        if (root[i])
        {
            for (int j = a; j < root.size(); j += 2 * i + 1) { root[j] = 0; }
        }
    }prime.push_back(2);
    for (int i = 1; i < root.size(); i++) { if (root[i] && jacobi(n, 2 * i + 1) == 1) { prime.push_back(2 * i + 1); } }
    root.clear();
    int m = 1 << (14 + std::min((int)d.num.size(), 3));
    for (int i = 1; i < prime.size(); i++)
    {
        int p = prime[i];
        logp.push_back(log(p));
        root.push_back(shanks_sqrt(n, p).num[0]);
    }
    d = sqroot(d * 14142 / (integer(m) * 10000));
    d.addsmall(31 - (d % 30).num[0]);
    q.p = n; q.init(1);
    smooth.assign(prime.size(), equa(n.num.size(), prime.size()));
    flaw.assign(prime.size() << 3, equa(n.num.size(), prime.size()));
    ok = 20;
    constexpr int z[] = { 6,4,2,4,2,4,6,2 }; int zi = 0;
    bound = getlog(d * m) * 2 - 13;
    while (ok > 0)
    {
        integer y0;
        do { d.addsmall(z[zi]); zi = (zi + 1) % 8; } while ((d.num[0] & 3) == 1 || (((y0 = power(n, (d + 1) / 4, d)) * y0 - n) % d).num.back());
        integer y1 = (n - y0 * y0) / d * inv(y0 * 2, d) % d;
        if (y1.sign == -1) { y1 = y1 + d; }
        integer a = d * d, b = y0 + y1 * d, c = (b * b - n) / a, e = b * inv(d, n) % n;
        std::vector<float>pos(m), neg(m);
        for (int i = 0; i < root.size(); i++)
        {
            int p = prime[i + 1], q = (a % p).num[0]; if (!q) { continue; }
            int h = p - 2, inva = h & 1 ? q : 1;
            while (h >>= 1)
            {
                q = q * q % p;
                if (h & 1) { inva = q * inva % p; }
            }
            h = (b % p).num[0];
            int x0 = (-h - root[i]) * inva % p;
            int x1 = (-h + root[i]) * inva % p;
            if (x0 < 0) { x0 += p; }if (x1 < 0) { x1 += p; }
            const float lp = logp[i];
            for (int j = x0; j < m; j += p) { pos[j] += lp; }
            for (int j = x1; j < m; j += p) { pos[j] += lp; }
            for (int j = p - x0; j < m; j += p) { neg[j] += lp; }
            for (int j = p - x1; j < m; j += p) { neg[j] += lp; }
        }
        for (int i = 0; i < m; i++)
        {
            if (pos[i] > bound)
            {
                g = insert_lambda(e + d * i, a * ((ll)i * i) + b * (i << 1) + c);
                if (check_lambda()) { return g; }
            }
            if (neg[i] > bound)
            {
                g = insert_lambda(e - d * i, a * ((ll)i * i) - b * (i << 1) + c);
                if (check_lambda()) { return g; }
            }
        }
    }
    std::cout << "factor"; exit(0);
}
integer ppow(const integer& x)
{
    integer p;
    if (x.num[0] % 2 == 0) { p = 2; }
    else if (x.num[0] % 5 == 0) { p = 5; }
    else
    {
        std::vector<int>b = tobinary(x);
        ll a = 2;
        mont q; q.p = x; q.init(1);
        do
        {
            p = gcd(q.out(q.pow_binary(q.in(a), b)) - a, x);
            if (a++ > 100 || (p.num.size() == 1 && p.num[0] == 1)) { return 0; }
        } while (!isprime(p));
    }
    std::vector<integer> l; l.push_back(p);
    while (x.absbigger(p = p * p, 1))
    {
        l.push_back(p);
    }
    p = x;
    do
    {
        integer r;
        p = p.divide(l.back(), r);
        if (r.num.back()) { return 0; }
        while (l.back().absbigger(p, 0))
        {
            if (l.size() == 1) { return p.num.size() == 1 && p.num[0] == 1 ? l[0] : 0; }
            l.pop_back();
        }
    } while (p.num.size() > 1 || p.num[0] > 1);
    return l[0];
}
integer euler(const integer& x, std::vector<integer>* l = 0, bool p = 1)
{
    if (x.num.size() == 1 && x.num[0] < 2) { return 1; }
    integer g = ppow(x);
    if (g.num.back()) { if (l) { l->push_back(g); }return x - x / g; }
    integer u = factor(x, p), v = x / u, n = u, m = v;
    integer k = g = gcd(u, v);
    while (g.num.size() > 1 || g.num[0] > 1)
    {
        u = u * g, v = v / g, g = gcd(u, v);
    }
    if (v.num.size() > 1 || v.num[0] > 1) { return euler(u, l, p) * euler(v, l, p); }
    u = m, v = n, g = k;
    while (g.num.size() > 1 || g.num[0] > 1)
    {
        u = u * g, v = v / g, g = gcd(u, v);
    }
    if (v.num.size() > 1 || v.num[0] > 1) { return euler(u, l, p) * euler(v, l, p); }
    return x / k * euler(k, l, p);
}
int mobius(const integer& x, bool p = 1)
{
    if (x.num.size() == 1 && x.num[0] == 1) { return 1; }
    integer g = ppow(x);
    if (g.num.back()) { return x.num == g.num ? -1 : 0; }
    integer u = factor(x, p), v = x / u;
    g = gcd(u, v);
    if (g.num.size() > 1 || g.num[0] > 1) { return 0; }
    int a = mobius(u, p);
    return a ? a * mobius(v, p) : 0;
}
int proot(const integer& p)
{
    if (p.num.size() == 1 && p.num[0] < 6) { return p.num[0] % 2 + 1; }
    std::vector<integer>k;
    integer x = p - 1;
    euler(x, &k);
    std::vector<integer>l;
    while (l.size() < k.size())
    {
        int t = -1;
        for (int i = 0; i < k.size(); i++)
        {
            if (k[i].sign > 0 && (t == -1 || k[t].absbigger(k[i], 0))) { t = i; }
        }
        l.push_back(k[t]); k[t].sign = -1;
    }k.clear();
    for (auto& e : l) { k.push_back(x / e); }l.clear();
    x = k[0];
    for (int i = 1; i < k.size(); i++)
    {
        l.push_back(k[i] - k[i - 1]);
    }
    mont q; q.p = p; q.init(1);
    for (int a = 2;; a++)
    {
        integer t = q.in(a), y = q.pow_binary(t, tobinary(x));
        int s = -1;
        while (y.num != q.f.num)
        {
            if (++s == l.size()) { return a; }
            y = q.out(y * q.pow_binary(t, tobinary(l[s])));
        }
    }
}
pairs fib(int n)
{
    pairs ans;
    if (n < 3)
    {
        ans.f = (n + 1) / 2, ans.s = (n + 2) / 2;
        return ans;
    }
    pairs  t = fib(n >> 1);
    integer a = t.f * t.f + t.s * t.s, b = t.f * (t.s + t.s - t.f);
    if (n & 1) { ans.f = a, ans.s = a + b; }
    else { ans.f = b, ans.s = a; }
    return ans;
}
#endif
