
//[Console]::OutputEncoding = [System.Text.Encoding]::UTF8
//g++ a.cpp -O2 -I D:\mingw64download\mingw64\include -L D:\mingw64download\mingw64\lib -lgmp -lgmpxx -o a.exe
#ifndef _INTEGERB_H_
#define _INTEGERB_H_
#include<iostream>
#include<vector>
#include<cmath>
#include<complex>
#define more 1//展示扩展部分算法
typedef long long ll;
typedef unsigned long long ull;
typedef std::complex<double> cd;
const int log2lenmax = 20, lenmax = 1 << log2lenmax, halflen = lenmax >> 1,e8=100000000;
constexpr int Blen=28,Base =1<<Blen,Bmask=Base-1,Flen=Blen>>1,Fmask=(1<<Flen)-1;
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
    static integer multiply(const int* a, int la, const int* b, int lb, int sign);
    static void div32(int*a,int bh,int bl,int&q)
    {
        ull p0=(ull)bl*q,p1=(ull)bh*q;
        a[0]-=p0&Bmask;
        bool carry=a[0]<0;
        a[0]+=carry<<Blen;
        a[1]-=(p0>>Blen)+(p1&Bmask)+carry;
        carry=a[1]<0;
        a[1]+=carry<<Blen;
        a[2]-=(p1>>Blen)+carry+(a[1]<0);
        a[1]+=(a[1]<0)<<Blen;
        while(a[2]<0)
        {
            q--;
            a[0]+=bl;
            carry=a[0]>=Base;
            a[0]-=carry<<Blen;
            a[1]+=bh+carry;
            carry=a[1]>=Base;
            a[1]-=carry<<Blen;
            a[2]+=carry;
        }
    }
    void quotient(const view& a, const view& b, integer& now)
    {
        int m =1,t = b.ptr[b.len - 1];
        while((t<<m)<Base){m++;}m-=1;
        int m_5=m+5;
        ull bb = ((ull)b.ptr[b.len - 1]<<m_5) + ( ((ull)b.ptr[b.len - 2]<<m)>>(Blen-5) );
        now.num.reserve(a.len + 1);
        now = a; now.num.push_back(0); int ln;
        integer prod; prod.num.resize(b.len + 2);
        if(a.len>b.len+30&&b.len>10)//朴素除法受进制影响极大,这里略微缓解，但是还是慢
        {
            int h=Blen-m,hmask=(1<<h)-1,bh=(b.ptr[b.len-1]<<m)|(b.ptr[b.len-2]>>h),
            bl=((b.ptr[b.len-2]&hmask)<<m)|(b.ptr[b.len-3]>>h);
            for (int tmp[4];(t = (ln = now.num.size())- b.len - 1) >10; now.num.pop_back())
            {
                int q = (((((ull)now.num.back()<<Blen) + now.num[ln - 2])<<m_5) + ( ((ull)now.num[ln - 3]<<m)>>(Blen-5) )) / bb;
                if (!q) { continue; }
                for(int i=3;i>-1;i--)
                {
                    tmp[i]=((now.num[ln-4+i]&hmask)<<m)|(now.num[ln-5+i]>>h);
                }
                div32(tmp+1,bh,bl,q);
                int q1=(((ull)tmp[2]<<Blen)+tmp[1])/bh;
                div32(tmp,bh,bl,q1);
                ull k=(ull)b.ptr[0] * q1;
                prod.num[0] = k&Bmask,k>>=Blen;
                for (int i = 1; i < b.len; i++)
                {
                    k = (ull)b.ptr[i] * q1+(ull)b.ptr[i-1]*q+ k;
                    prod.num[i] = k&Bmask,k>>=Blen;
                }
                k=(ull)b.ptr[b.len-1]*q+k;
                prod.num[b.len] = k&Bmask,k>>=Blen;
                prod.num[b.len+1] = k;
                view pview(prod.num.data(), 0, b.len+1, 1);
                int s = abssub(now.num.data() + t-1, b.len + 2, pview.ptr, pview.len);
                int cnt = 0;
                while (s < 0)
                {
                    s = -abssub(now.num.data() + t-1, b.len +2, b.ptr, b.len);
                    if (cnt++ > 2) {std::cout << "quotient1"; exit(0); }
                }
                ull q2=(((ull)q<<Blen)|q1)-cnt;
                num[t]=q2>>Blen;
                num[t-1]=q2&Bmask;
                now.num.pop_back();
            }
        }
        for (; (t = (ln = now.num.size()) - b.len - 1) > -1; now.num.pop_back())
        {
            int q = (((((ull)now.num.back() << Blen) + now.num[ln - 2]) << m_5) + (((ull)now.num[ln - 3] << m) >> (Blen - 5))) / bb;
            if (!q) { continue; }
            ull k = 0;
            for (int i = 0; i < b.len; i++)
            {
                k = (ull)b.ptr[i] * q + k;
                prod.num[i] = k & Bmask, k >>= Blen;
            }prod.num[b.len] = k;
            view pview(prod.num.data(), 0, b.len, 1);
            int s = abssub(now.num.data() + t, b.len + 1, pview.ptr, pview.len);
            int cnt = 0;
            while (s < 0)
            {
                s = -abssub(now.num.data() + t, b.len + 1, b.ptr, b.len);
                if (cnt++> 2) { std::cout << "quotient2 ";exit(0); }
            }
            num[t] = q - cnt;
        }
        while (now.num.back() == 0 && now.num.size() > 1) { now.num.pop_back(); }
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
public:
    static int abssub(int* a, int la, const int* b, int lb);
    static integer div_native(const view& a, const view& b, integer& r)
    {
        if (!b.ptr[b.len - 1]) { std::cout << "div0"; exit(0); }
        integer c;
        c.num.resize(std::max(a.len - b.len + 1, 1));
        if (b.len == 1)//now*base不溢出
        {
            ull n = b.ptr[0], now = 0;
            for (int i = a.len - 1; i > -1; i--)
            {
                now = (now<<Blen) + a.ptr[i];
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
            num.push_back(n&Bmask);
            n>>=Blen;
        } while (n);
    }
    void get()
    {
        std::string s;
        std::cin >> s;
        *this = integer(s);
    }
    integer(const std::string& s)
    {
        int i = s.size() - 8;
        for (; i > 0; i -= 8)
        {
            int a = s[i] - 48;
            for (int j = 1; j < 8; j++)
            {
                a *= 10;
                a += s[i + j] - 48;
            }
            num.push_back(a);
        }
        i += 8;
        int a = 0, j = 0;
        if (s[0] == '-') { j = 1; sign = -1; }
        for (; j < i; j++)
        {
            a *= 10;
            a += s[j] - 48;
        }
        if (a || num.empty()) { num.push_back(a); }
        int n=1;while(n<num.size()){n<<=1;}
        num.resize(n);
        integer h;h.num.push_back(e8);
        for(int len=1;len<n;)
        {
            int len2=len<<1;
            for(int i=0;i<n;i+=len2)
            {
                view high(num.data()+i,len,len2-1,1);
                view low(num.data()+i,0,len-1,1);
                integer tmp=karamul(high,h);
                tmp=addorsub(tmp.num.data(),tmp.num.size(),1,low.ptr,low.len,1,1);
                int j=0;
                for(;j<tmp.num.size();j++)
                {
                    num[i+j]=tmp.num[j];
                }
                for(;j<len2;j++){num[i+j]=0;}
            }
            len=len2;
            if(len<n){h=h*h;}
            else{break;}
        }
        while(num.back()==0&&num.size()>1){num.pop_back();}
    }
    ll toll()const
    {
        int ln = num.size();
        if (ln > 2) { std::cout << "toll"; exit(0); }
        ll a = num.back();
        if (ln > 1) { a = (a<<Blen)+ num[0]; }
        return sign * a;
    }
    std::string tostring(int m=0)const
    {
        if(m<1){m=num.size()<<2;}//全精度
        std::string s;
        if (sign != 1) { s.push_back('-'); }
        if(num.back()){
            integer E8;E8.num.push_back(e8);
            int l = num.size(),n=0;integer tmp;
            if(l>m)
            {
                n=(l-m)*(log(Base)/log(e8))-1;
            }
            if(n<1)
            {
                n=0;
                tmp=*this;
            }
            else{tmp=*this/shiftpow(E8,n,num.size(),l);}
            if(tmp.absbigger(E8,1))
            {
                std::vector<integer>p;p.push_back(E8);
                int t=2;
                while(tmp.absbigger(p.back(),0)&&p.back().num.size()*2-1<=tmp.num.size())
                {
                    t<<=1;p.push_back(p.back()*p.back());
                }
                if(p.back().absbigger(tmp,0)){p.pop_back();t/=2;}
                tmp.num.resize(t);
                for(int len=t;len>1;)
                {
                    int len2=len>>1;
                    for(int i=0;i<t;i+=len2)
                    {
                        view a(tmp.num.data()+i,0,len-1,1);
                        integer r,q=divide(a,p.back(),r);
                        int j=0;
                        for(;j<r.num.size();j++)
                        {
                            tmp.num[i+j]=r.num[j];
                        }
                        for(;j<len2;j++){tmp.num[i+j]=0;}
                        j=0,i+=len2;
                        for(;j<q.num.size();j++)
                        {
                            tmp.num[i+j]=q.num[j];
                        }
                        for(;j<len2;j++){tmp.num[i+j]=0;}
                    }
                    len=len2;
                    if(len>1){p.pop_back();}
                    else{break;}
                }
                while(tmp.num.back()==0){tmp.num.pop_back();}
            }
            s.reserve(tmp.num.size()*8+20);
            int cnt = 0; int a[8];
            a[7] = tmp.num.back();
            for (int j = 7; j > 0; j--)
            {
                a[j - 1] = a[j] / 10;
            }
            for (; cnt < 7 && a[cnt] == 0; cnt++);
            if (n) { 
                s.push_back('0' + (a[cnt++] % 10));
                s.push_back('.'); 
            }
            for (int j = cnt; j < 8; j++)
            {
                s.push_back('0' + (a[j] % 10));
            }
            for (int i = 1; i < tmp.num.size(); i++)
            {
                a[7] = tmp.num[tmp.num.size() - 1 - i];
                for (int j = 7; j > 0; j--)
                {
                    a[j - 1] = a[j] / 10;
                }
                for (int j = 0; j < 8; j++)
                {
                    s.push_back('0' + (a[j] % 10));
                }
            }
            if (n)
            {
                s.push_back('e');
                s += std::to_string((tmp.num.size() + n) * 8 - cnt);
            }
        }
        else{
            s.push_back('0');
        }
        return s;
    }
    void  print(bool stdendl=1,int m=5)const
    {
        std::cout<<tostring(m);
        if (stdendl) { std::cout<<'\n'; }
        else { std::cout<<' '; }
    }
    bool absbigger(const view& that, bool equal) const
    {
        int n = num.size();
        if (n == that.len)
        {
            int i = n - 1;
            for (; i > -1 && num[i] == that.ptr[i]; i--);
            if (i == -1) { return equal; }
            return num[i] > that.ptr[i];
        }
        return n > that.len;
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
    static integer shiftmul(const view& a, const view& b, int need, int& offset)
    {
        if (need < 1) { std::cout << "shiftmul"; exit(0); }
        int b1 = std::max(0, a.len - need);
        int b2 = std::max(0, b.len - need);
        offset += b1 + b2;
        view aview(a, b1), bview(b, b2);
        return karamul(aview, bview);
    }
    static integer fftmul(const view& a, const view& b);
    static double* b14rrii(const view& a, int n);
    integer operator*(const integer& that)const
    {
        return karamul(*this, that);
    }
    static void reciprocal(const view& n, integer& xt, int l)//xt=base^t/*this
    {
        l++;
        int k = 1;
        int t = n.len + k;
        int b;
        integer tmp;
        xt = div_native(integer(1).shift(3), view(n.ptr, n.len - 2, n.len - 1, n.sign), tmp);
        std::vector<int>xtsize;
        int k1 = l - xt.num.size() + 1; k1 += ((k1 & 1) == 0);
        while (k1 > 5) { k1 = k1 / 2 + 1; k1 += ((k1 & 1) == 0); xtsize.push_back(k1); }
        while (t < l) {
            if (t + k >= l) { k = l - t; }
            else if (!xtsize.empty() && xtsize.back() > xt.num.size() && xtsize.back() - xt.num.size() < k) { k = xtsize.back() - xt.num.size(); xtsize.pop_back(); }
            int need = xt.num.size() + k + 2;
            b = 0;
            tmp = shiftmul(xt, xt, need, b);
            tmp = shiftmul(n, tmp, need, b);
            xt = (xt + xt).shift(k);
            //std::cout<<xt.num.back()<<" "<<xt.num.size()<<" "<< -(b + k - t)<<"\n";
            view tmpview(tmp, -(b + k - t));
            abssub(xt.num.data(), xt.num.size(), tmpview.ptr, tmpview.len);
            while (xt.num.back() == 0 && xt.num.size() > 1) { xt.num.pop_back(); }
            t += k;
            k = xt.num.size() - 1;
            while (!xtsize.empty() && xt.num.size() >= xtsize.back()) { xtsize.pop_back(); }
        }
        xt = xt.shift(-1);
        xt.addsmall(-n.sign);
    }
    integer fsqrt()const//xt=sqrt(*this/base^(2f-2k))
    {
        int ns = num.size();
        if (ns< 3)
        {
            ll a = toll(), b = sqrt(a);
            if(ns<2){return b;}
            b=(b+1+a/(b+1))/2;
            while(b * b > a) { b--; }
            return b;
        }
        int k = 2;
        int f=ns/2+(ns%2!=0);
        const bool useyk = f > 100;
        int need, b;
        integer xt;
        double guess= (double)((ll)num.back() * Base + num[ns - 2])*Base+num[ns-3];
        if (ns % 2)
        {
            xt = sqrt(guess)+1;
        }
        else
        {
            xt = sqrt(guess*Base+num[ns-4])+1000;
            integer q;need=k+1;
            do
            {
                b=0;
                integer tmp=shiftmul(xt,xt,need,b).shift(b);
                view nview(*this,2*(f-k));
                abssub(tmp.num.data(),tmp.num.size(),nview.ptr,nview.len);
                while(tmp.num.size()>1&&tmp.num.back()==0){tmp.num.pop_back();}
                q=tmp/xt;q.div2pow(1);
            } while(q.sign>0&&q.num.back()&&(xt=xt-q).num.back());
        }
        int p=1;
        integer yt; int l,yk=0;
        std::vector<int>kplan;
        int k1 = f; k1 += ((k1 & 1) == 0);
        while (k1 > 5) { k1 = k1 / 2 + 1; k1 += ((k1 & 1) == 0); kplan.push_back(k1); }
        while (k<f)
        {
            //std::cout<<k<<" "<<f<<"\n";
            if (p+ k >f) { p=f-k; }
            else if (!kplan.empty() && kplan.back() > k&& kplan.back() - k <p) 
            { p = kplan.back() - k; kplan.pop_back(); }
            need = xt.num.size() + p + 2, b = 0;
            integer tmp = shiftmul(xt, xt, need, b).shift(p+ b);
            view nview(*this, 2*(f-k)-p);
            int check = abssub(tmp.num.data(), tmp.num.size(), nview.ptr, nview.len);
            if (check < 0) { std::cout << "fsqrt"; exit(0); }
            while (tmp.num.size() > 1 && tmp.num.back() == 0) { tmp.num.pop_back(); }
            integer q;bool flag=useyk;
            if (useyk)
            {
                if (yt.num.empty())
                {
                    l = 2*xt.num.size();
                    yt = integer(1).shift(l) / xt;
                }
                else
                {
                    yk = (int)xt.num.size() - (int)yt.num.size();
                    need = yt.num.size() + yk + 2, b = 0;
                    integer ytmp = shiftmul(yt, yt, need, b);
                    ytmp = shiftmul(xt, ytmp, need, b);
                    view ytmpview(ytmp, -(yk - l + b));
                    yt = (yt + yt).shift(yk);
                    abssub(yt.num.data(), yt.num.size(), ytmpview.ptr, ytmpview.len);
                    while (yt.num.back() == 0 && yt.num.size() > 1) { yt.num.pop_back(); }
                    l += yk;
                }
                need = tmp.num.size() - xt.num.size() + 1, b = 0;
                if(need>7)
                {
                    q = shiftmul(tmp, yt,need, b);
                    view qview(q, l - b);
                    integer prod = karamul(qview, xt);
                    tmp.sign = abssub(tmp.num.data(), tmp.num.size(), prod.num.data(), prod.num.size());
                    while (tmp.num.back() == 0 && tmp.num.size() > 1) { tmp.num.pop_back(); }
                    //std::cout<<tmp.num.size()<<" "<<xt.num.size()<<"\n";
                    integer r;
                    tmp = divide(tmp,xt,r);
                    q = addorsub(qview.ptr, qview.len, qview.sign, tmp.num.data(), tmp.num.size(), tmp.sign, 1);
                    if(tmp.sign==-1&&r.num.back()){q.addsmall(-1);}
                }else{flag=0;}
            }
            if(!flag)
            {
                q = tmp / xt;
            }
            q.div2pow(1);
            xt = xt.shift(p);
            abssub(xt.num.data(), xt.num.size(), q.num.data(), q.num.size());
            while (xt.num.size() > 1 && xt.num.back() == 0) { xt.num.pop_back(); }
            l +=p;
            k+=p;
            p = xt.num.size() - 1;
            while (!kplan.empty() && k >= kplan.back()) { kplan.pop_back(); }
        }
        xt.addsmall(-1);
        integer prod = xt * xt;
        if (prod.absbigger(*this, 0)){xt.addsmall(-1);}
        return xt;
    }
    integer root(int m)const
    {  
        if (!num.back() || m == 1) { return *this; }
        if (sign < 0 || m < 1) { std::cout << "root"; exit(0); }
        if (m == 2) { return fsqrt(); }
        int ns = num.size();
        int f=ns/m+(ns%m!=0);
        const bool useyk =f>100;
        int k=1+(f>1),need=k+1, b,b1;
        view h(*this,m*(f-k));
        double guess=h.len== 1 ? log(h.ptr[0]) : log(((ll)h.ptr[h.len-1]<<Blen) + h.ptr[h.len- 2]) + (h.len- 2) * log(Base);
        integer target=exp(guess/m);
        integer xt=integer(1).shift(k),q=integer(1).shift(k-1);
        int gap=f>1?Base/10:1;
        ll rate=(ll)1<<50;
        if(q.absbigger(target,1)){target=q+1;rate>>=10;}
        else if(target.absbigger(xt,1)){target=xt-1;rate>>=10;}
        while(xt.absbigger(q+gap,0))
        {
            integer mid=(xt+q)/2;
            if(rate<2){rate=2;}
            mid=(target*(rate-1)+mid)/rate;
            b=0;
            integer tmp=shiftpow(mid,m,need,b);
            int leftshift=m*(f-k)+b;
            if(leftshift<ns&&!tmp.absbigger(view(*this,leftshift),0))
            {
                q=mid;
                if(q.absbigger(target,1)){target=q+1;rate>>=10;}
            }
            else
            {
                xt=mid;
                if(target.absbigger(xt,1)){target=xt-1;rate>>=10;}
            }
        }
if(f>1)
{
        do
        {
            b=0;
            integer tmp=shiftpow(xt,m-1,need,b);
            b1=0;
            integer tmp1=shiftmul(xt,tmp,need,b1).shift(b1);
            int leftshift=m*(f-k)+b;
            if(ns>leftshift){
            view nview(*this,leftshift);
            abssub(tmp1.num.data(),tmp1.num.size(),nview.ptr,nview.len);
            while(tmp1.num.size()>1&&tmp1.num.back()==0){tmp1.num.pop_back();}
            }
            q=tmp1/(tmp*m);
        } while(q.sign>0&&q.num.back()&&(xt=xt-q).num.back());
        int p=k-1;
        integer yt; int l,yk=0;
        std::vector<int>kplan;
        int k1 = f; k1 += ((k1 & 1) == 1);
        while (k1 > 5) { k1 = k1 / 2+1; k1 += ((k1 & 1) == 1); kplan.push_back(k1); }
        int cnt=0;//必须保护第一轮
        while (k<f||cnt==1)
        {
            //std::cout<<k<<"\n";
            cnt++;
            if (p+ k >f) { p=f-k; }
            else if (!kplan.empty() && kplan.back() > k&& kplan.back() - k <p) { p = kplan.back() - k; kplan.pop_back(); }
            need = xt.num.size() + p + 2, b = 0,b1=0;
            integer tmp=shiftpow(xt,m-1,need,b);
            integer tmp1=shiftmul(xt,tmp,need,b1).shift(b1+p);
            int leftshift=m*(f-k)+b-p;
            if(ns>leftshift){
            view nview(*this,leftshift);
            int check = abssub(tmp1.num.data(), tmp1.num.size(), nview.ptr, nview.len);
            if (check < 0) { std::cout << "root"; exit(0); }
            while (tmp1.num.size() > 1 && tmp1.num.back() == 0) { tmp1.num.pop_back(); }
            }
            bool flag=useyk;
            if (useyk)
            {
                int diff=tmp.num.size()-xt.num.size();
                view tmpview(tmp,diff);//与xt.size对齐
                if (yt.num.empty())
                {
                    l = 2*tmpview.len;
                    integer r;
                    yt = divide(integer(1).shift(l),tmpview,r);
                }
                else
                {
                    yk = tmpview.len- (int)yt.num.size();
                    need = yt.num.size() + yk + 2, b = 0;
                    integer ytmp = shiftmul(yt, yt, need, b);
                    ytmp = shiftmul(tmpview, ytmp, need, b);
                    view ytmpview(ytmp, -(yk - l + b));
                    yt = (yt + yt).shift(yk);
                    abssub(yt.num.data(), yt.num.size(), ytmpview.ptr, ytmpview.len);
                    while (yt.num.back() == 0 && yt.num.size() > 1) { yt.num.pop_back(); }
                    l += yk;
                }
                need = tmp1.num.size() -tmp.num.size() + 1, b = 0;
                if(need>7)
                {
                    q = shiftmul(tmp1, yt,need, b);
                    view qview(q, l+diff- b);
                    integer prod = karamul(qview,tmp);
                    tmp1.sign= abssub(tmp1.num.data(), tmp1.num.size(), prod.num.data(), prod.num.size());
                    while (tmp1.num.back() == 0 && tmp1.num.size() > 1) { tmp1.num.pop_back(); }
                    //std::cout<<tmp.num.size()<<" "<<tmp1.num.size()<<"\n";
                    integer r;
                    tmp = divide(tmp1,tmp,r);
                    q = addorsub(qview.ptr, qview.len, qview.sign, tmp.num.data(), tmp.num.size(), tmp.sign, 1);
                    if(tmp.sign==-1&&r.num.back()){q.addsmall(-1);}
                    q = q /m;
                }else{flag=0;}
            }
            if(!flag)
            {
                q = (tmp1/tmp) /m;
            }
            xt = xt.shift(p);
            abssub(xt.num.data(), xt.num.size(), q.num.data(), q.num.size());
            while (xt.num.size() > 1 && xt.num.back() == 0) { xt.num.pop_back(); }
            l +=p;
            k+=p;
            p = xt.num.size() -2;//必须更加保守
            while (!kplan.empty() && k >= kplan.back()) { kplan.pop_back(); }
        }
    }
        xt.addsmall(-1);
        need=num.size()<<1;//全精度
        if(shiftpow(xt,m,need,b).absbigger(*this,0)){xt.addsmall(-1);}
        return xt;
    }
    static integer shiftpow(const integer& x, int n, int need, int& b)
    {
        integer a(x), ans = n & 1 ? x : 1;
        int b1 = 0;
        while (n >>= 1)
        {
            b1 <<= 1;
            a = shiftmul(a, a, need, b1);
            if (n & 1) { ans = shiftmul(ans, a, need, b); b += b1; }
        }
        return ans;
    }
    static integer div_newton(const view& a, const view& b, integer& r)
    {
        r.num.clear();
        int la = a.len, lb = b.len;
        int n = std::min(24, la / lb);
        int l = (la - lb) / n + 3;
        while (l >= 512) { l >>= 1; }n += l < 285 && l > 255;
        l = (la - lb) / n + lb + 1;
        integer xt; reciprocal(b, xt, l);
        integer res(a);
        integer q; q.num.assign(la - lb + 1, 0);
        for (int i = 0; i < n && res.num.size() > lb + 4; i++) {
            int lr = res.num.size(), start = lr - l;
            if (start < 0) { start = 0; }
            view cutview(res.num.data(), start, lr - 1, res.sign);
            int need = cutview.len - lb + 2, b_shift = 0;
            integer p = shiftmul(xt, cutview, need, b_shift);
            view pview(p, l - b_shift);
            integer now = karamul(b, pview);
            if (abssub(res.num.data() + start, std::min(l, lr), now.num.data(), now.num.size()) < 0) { std::cout << "div_newton1"; exit(0); }
            q.shiftadd(pview, start);
            while (!res.num.back() && res.num.size() > 1) { res.num.pop_back(); }
        }
        if (res.num.size() > lb + 4) { std::cout << "div_newton2"; exit(0); }
        q.shiftadd(div_native(res, b, r), 0);
        while (q.num.size() > 1 && q.num.back() == 0) { q.num.pop_back(); }
        q.sign = a.sign * b.sign;
        r.sign = a.sign;
        return q;
    }
    static integer divide(const view& a, const view& b, integer& r)
    {
        int la = a.len, lb = b.len, c = la - lb;
        if (la>200&&lb>50&&c >20&&(c>100||lb * c >10000))
        {
            return div_newton(a, b, r);
        }
        return div_native(a, b, r);
    }
    integer operator/(const integer& that)const
    {
        integer r;
        return divide(*this, that, r);
    }
    integer operator%(const integer& that)const//符号同被除数
    {
        integer r;
        divide(*this, that, r);
        return r;
    }
    integer mod_positive(const integer&b)const//返回非负
    {
        integer r;
        divide(*this,b, r);
        if(r.sign<0&&r.num.back())
        {
            r.num.resize(b.num.size());//在divide已经reserve足够空间
            abssub(r.num.data(),r.num.size(),b.num.data(),b.num.size());
            while(r.num.back()==0&&r.num.size()>1){r.num.pop_back();}
            r.sign=1;
        }
        return r;
    }
    int ctz()const
    {
        int i=0;
        for(;i<num.size()-1&&num[i]==0;i++);
        int h=num[i];i*=Blen;
        while((h&1)==0){h>>=1;i++;}
        return i;
    }
    void mod2pow(int n)
    {
        int end=n/Blen;n%=Blen;
        if(end>=num.size()){return;}
        num[end]&=(1<<n)-1;
        num.resize(end+1);
    }
    void div2pow(int n)
    {
        int start=n/Blen,i=0;n%=Blen;
        int m=Blen-n,mask=(1<<n)-1;
        for(;start+1<num.size();i++)
        {
            int low=num[start]>>n,high=num[++start]&mask;
            num[i]=low|(high<<m);
        }
        if(start>=num.size()){num[i]=0;}
        else
        {
            int h=num[start]>>n;
            if(h){num[i]=h;}
            else{num[i]=0;i-=i>0;}
        }
        num.resize(i+1);
    }
    void mul2pow(int n)
    {
        if(!num.back()){return;}
        int start=n/Blen;n%=Blen;
        int m=Blen-n,mask=(1<<m)-1,ns=num.size();
        num.resize(start+ns+1);
        for(int i=ns;i>0;i--)
        {
            int high=num[i]&mask,low=num[i-1]>>m;
            num[start+i]=low|(high<<n);
        }
        num[start]=(num[0]&mask)<<n;
        for(int i=0;i<start;i++){num[i]=0;}
        if(num.back()==0){num.pop_back();}
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
            nextlevel.reserve(level.size());
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
            a[i] += k<<Blen;
        }
    }
    else
    {
        int i = 0; s = 1;
        for (; i < lb && i < la; i++)
        {
            a[i] -= k;
            k = a[i] < 0;
            a[i] += k<<Blen;
        }
        for (; i < la && k; i++)
        {
            a[i] -= k;
            k = a[i] < 0;
            a[i] += k<<Blen;
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
            result.num[i] -= k<<Blen;
        }
        const int lr = result.num.size();
        for (j = lb; j < lr && k; j++)
        {
            result.num[j] += k;
            k = result.num[j] >= Base;
            result.num[j] -= k<<Blen;
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
    if (!lb||(la + a[0] < 2) || (lb + b[0] < 2)) { return 0; }
    int n = la - 1, m = lb - 1, l = n + m + 1;
    integer c;
    c.sign = sign;
    c.num.resize(l + 1);
    ull k = 0;
    int i;
    for (i = 0; i < m; i++)
    {
        int j = i;
        for (; j > 3; j -= 4)
        {
            int t = i - j;
            k += (ll)a[t] * b[j]
                + (ll)a[t + 1] * b[j - 1]
                + (ll)a[t + 2] * b[j - 2]
                + (ll)a[t + 3] * b[j - 3];
        }
        for (; j > -1; j--)
        {
            k += (ll)a[i - j] * b[j];
        }
        c.num[i] = k&Bmask;
        k>>=Blen;
    }
    for (; i < n; i++)
    {
        int j = m;
        for (; j > 3; j -= 4)
        {
            int t = i - j;
            k += (ll)a[t] * b[j]
                + (ll)a[t + 1] * b[j - 1]
                + (ll)a[t + 2] * b[j - 2]
                + (ll)a[t + 3] * b[j - 3];
        }
        for (; j > -1; j--)
        {
            k += (ll)a[i - j] * b[j];
        }
        c.num[i] = k&Bmask;
        k>>=Blen;
    }
    for (; i < l; i++)
    {
        int j = m;
        for (; j > i - n + 3; j -= 4)
        {
            int t = i - j;
            k += (ll)a[t] * b[j]
                + (ll)a[t + 1] * b[j - 1]
                + (ll)a[t + 2] * b[j - 2]
                + (ll)a[t + 3] * b[j - 3];
        }
        for (; j > i - n - 1; j--)
        {
            k += (ll)a[i - j] * b[j];
        }
        c.num[i] = k&Bmask;
        k>>=Blen;
    }
    if(k){c.num[l] = k;}
    else{ c.num.pop_back(); }
    return c;
}
integer integer::karamul(view a, view b)
{
    if (a.len < b.len) { std::swap(a, b); }
    if (b.len < 80 || (b.len < 160 && a.len * b.len < 30000)) {return multiply(a.ptr, a.len, b.ptr, b.len, a.sign * b.sign);}
    if (a.len + b.len <= halflen) {return fftmul(a, b);}
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
        ans.num = c0.num;
        ans.num.resize(a.len + b.len);
        integer tmp = karamul(addorsub(a0.ptr, a0.len, 1, a1.ptr, a1.len, 1, 1), addorsub(b0.ptr, b0.len, 1, b1.ptr, b1.len, 1, 1));
        abssub(tmp.num.data(), tmp.num.size(), c0.num.data(), c0.num.size());
        abssub(tmp.num.data(), tmp.num.size(), c1.num.data(), c1.num.size());
        while (tmp.num.back() == 0 && tmp.num.size() > 1) { tmp.num.pop_back(); }
        ans.shiftadd(tmp, n);
        ans.shiftadd(c1, 2 * n);
    }
    while (ans.num.size() > 1 && ans.num.back() == 0) ans.num.pop_back();
    ans.sign = a.sign * b.sign;
    return ans;
}
double* integer::b14rrii(const view& a, int n)
{
    double* y = new double[n]();
    int i = 0;
    for (; i + 1 < a.len; i += 2)
    {
        y[2 * i] = a.ptr[i]&Fmask;
        y[2 * i + 2] = a.ptr[i]>>Flen;
        y[2 * i + 1] = a.ptr[i + 1]&Fmask;
        y[2 * i + 3] = a.ptr[i + 1]>>Flen;
    }
    if (i != a.len)
    {
        y[2 * i] = a.ptr[i]&Fmask;
        y[2 * i + 2] = a.ptr[i]>>Flen;
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
    double* x0 = b14rrii(a, n), * y0;
    c2* xx = reinterpret_cast<c2*>(&x0[0]), * yy;
    fft(xx, n4);
    const bool same = a.ptr == b.ptr && a.len == b.len;
    if (same) { y0 = x0, yy = xx; }
    else { y0 = b14rrii(b, n), yy = reinterpret_cast<c2*>(&y0[0]); fft(yy, n4); }
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
    c.num.reserve(m = la + lb + 1); m /= 2;
    c.sign = a.sign * b.sign;
    ll k = 0;
    for (int i = 0; i < m; i++)
    {
        ll t = (ll)(xx[i].r.a + 0.5) + k;
        ll pre = t&Fmask;
        k = t>>Flen;
        t = (ll)(xx[i].i.a + 0.5) + k;
        k = t>>Flen;
        c.num.push_back(((t&Fmask)<<Flen) + pre);
        t = (ll)(xx[i].r.b + 0.5) + k;
        pre = t&Fmask;
        k = t>>Flen;
        t = (ll)(xx[i].i.b + 0.5) + k;
        k = t>>Flen;
        c.num.push_back(((t&Fmask)<<Flen) + pre);
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
    }
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
    b = cnt(e, 2) - cnt(s, 2) - cnt(e_s, 2);
    for (int i = 0, j = 0; i < u.size() - 2; j = u[i], i += 2)
    {
        l.push_back(power(mul_lambda(j, u[i]), u[i + 1]));
    }
    integer::mul_core(l, s);
    integer r = (l.back() * power(2, b % Blen));
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
    std::vector<integer> level;level.reserve(s/2+10);
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
    ll xh=x.num.back(),yh=y.num.back();
    if(x.num.size()>1)
    {
        xh=(xh<<Blen)+x.num[x.num.size() - 2];
        if(x.num.size()==y.num.size())
        {
            yh = (yh<<Blen) + y.num[y.num.size() - 2];
        }
    }
    ll A = 1, B = 0, C = 0, D = 1, q, u, v,check;
    while ((u = yh + B) && (v = yh + D) && (check=xh+C-(q = (xh + A) / u)*v)>-1&&check<v)
    {
        u = A, A = B, B = u - B * q;
        u = C, C = D, D = u - C * q;
        u = xh, xh = yh, yh = u - yh * q;
    }
    if (B)
    {
        int i = 0; ll kx = 0, ky = 0;
        for (; i < y.num.size(); i++) {
            int xi = x.num[i], yi = y.num[i];
            kx += A * xi + C * yi;
            ky += B * xi + D * yi;
            x.num[i] = kx&Bmask; kx>>=Blen;//kx负数也行
            y.num[i] = ky&Bmask; ky>>=Blen;
        }
        for (; i < x.num.size(); i++) {
            kx += A * x.num[i];
            x.num[i] = kx&Bmask; kx>>=Blen;
        }
        while (x.num.size() > 1 && x.num.back() == 0) x.num.pop_back();
        while (y.num.size() > 1 && y.num.back() == 0) y.num.pop_back();
        auto shift2 = [](integer& a, integer& b, ll A, ll B, ll C, ll D) {
            int i = 0; ll ka = 0, kb = 0;
            if (!b.num.back()) b.sign = -a.sign;
            for (; i < a.num.size(); i++) {
                int ai = a.num[i], bi = b.num[i];
                ka += llabs(A * ai - C * bi);
                kb += llabs(B * ai - D * bi);
                a.num[i] = ka&Bmask; ka>>=Blen;
                b.num[i] = kb&Bmask; kb>>=Blen;
            }
            for (; i < b.num.size(); i++) {
                ka += llabs(C * b.num[i]);
                kb += llabs(D * b.num[i]);
                a.num.push_back(ka&Bmask); ka>>=Blen;
                b.num[i] = kb&Bmask; kb>>=Blen;
            }
            while (ka) {  a.num.push_back(ka&Bmask); ka>>=Blen; }
            while (kb) {  b.num.push_back(kb&Bmask); kb>>=Blen; }
            b.sign *= (D > 0 ? 1 : -1);
            a.sign = -b.sign;
            };
        if (a) shift2(*a, *b, A, B, C, D), shift2(*c, *d, A, B, C, D);
    }
    return B;
}
void hgcd(const integer::view& x, const integer::view& y, integer& a, integer& b, integer& c, integer& d)//多项式版本粗暴映射
{
    int m = x.len / 2;
    if (m < 100)
    {
        integer xn(x), yn(y);
        while (yn.num.back() && yn.num.size() >= m)
        {
            int l = yn.num.size();
            if (!(l > m && xn.num.size() < l + 2 && lehmer(xn, yn, &a, &b, &c, &d)))
            {
                integer r, q = integer::divide(xn, yn, r);
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
    A = 0, B = C = 1, D = integer::divide(xn, yn, r); D.sign = -D.sign;
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
        while (y.num.size() > 3100 && 10 * y.num.size() > 9 * x.num.size())
        {
            integer A(1), B(0), C(0), D(1);
            hgcd(x, y, A, B, C, D);
            gcdshift1(x, y, A, B, C, D);
        }
        if (y.num.back())
        {
            if (!(x.num.size()< y.num.size() + 2 && lehmer(x, y, 0, 0, 0, 0)))
            {
                integer r; integer::divide(x, y, r);
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
    const bool big=y.num.size()>900;
    while (y.num.size() > 900)
    {
        while (10 * y.num.size() > 9 * x.num.size())
        {
            integer A(1), B(0), C(0), D(1);
            hgcd(x, y, A, B, C, D);
            gcdshift1(x, y, A, B, C, D);
            gcdshift2(a, b, c, d, A, B, C, D);
        }
        if (y.num.size() > 900)
        {
            integer r, q = integer::divide(x, y, r);
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
            if (!(x.num.size()< y.num.size() + 2 && lehmer(x, y, &A, &B, &C, &D)))
            {
                integer r, q = integer::divide(x, y, r);
                gcdshift(A, B, q);
                gcdshift(C, D, q);
                std::swap(x, y); std::swap(y, r);
            }
        }
        if(big){a=a*A+b*C,c=c*A+d*C;}
        else{std::swap(a,A),std::swap(c,C);}
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
    return t.mod_positive(p);
}
class mont
{
private:
    integer a, b;
    int m;
    integer::view fastmod(const integer& x)
    {
        return  integer::view(x.num.data(), 0, std::min((int)x.num.size() - 1, m - 1), x.sign);
    }
    integer::view fastdiv(const integer& x)
    {
        return integer::view(x.num.data(),m,x.num.size()-1,1);
    }
public:
    integer f, p;
    mont() {}
    void init()
    {
        m = p.num.size();
        integer n= integer(1).shift(m);
        a = n - inv(p, n);
        f = n % p;
        b = f * f % p;
    }
    integer out(const integer& x)
    {
        integer tmp1=integer::karamul(fastmod(x),a);
        integer tmp2=x+integer::karamul(p,fastmod(tmp1));
        if(tmp2.num.size()<=m){return 0;}
        integer::view t = fastdiv(tmp2);
        if (!p.absbigger(t, 0)) 
        { 
            integer::abssub(tmp2.num.data()+m,t.len,p.num.data(),p.num.size());
            while(t.len>1&&t.ptr[t.len-1]==0){t.len--;} 
        }
        return t;
    }
    integer in(const integer& x)
    {
        return out((x * b).mod_positive(p));
    }
    integer add(const integer& u, const integer& v)
    {
        integer ans(u + v);
        if (ans.absbigger(p, 1)) 
        { 
            integer::abssub(ans.num.data(),ans.num.size(),p.num.data(),p.num.size());
            while(ans.num.size()>1&&ans.num.back()==0){ans.num.pop_back();} 
        }
        return ans;
    }
    integer sub(const integer& u, const integer& v)
    {
        integer ans;
        if(u.absbigger(v,1)){ans=u;}
        else{ans=u+p;}
        integer::abssub(ans.num.data(),ans.num.size(),v.num.data(),v.num.size());
        while(ans.num.size()>1&&ans.num.back()==0){ans.num.pop_back();} 
        return ans;
    }
    integer pow_binary(const integer& x, const integer&y)
    {
        const std::vector<int>&l=y.num;
        if (l.size() == 1 && l[0] == 0) { return f; }
        integer u(x), v;
        bool k = 1;
        const int ls = l.size() - 1;
        for (int i = 0; i <= ls; i++)
        {
            ll m = l[i];
            for (int j = 0; j < 28; j++)
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
    integer pow_exponent(const integer& u, const std::vector<int>& l, int w)
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
    int a1,m;
    integer pow_exponent_m(const integer& u, const std::vector<int>& l,int w)
    {
        if (l.empty()) { return 1; }
        std::vector<integer>b; b.push_back(u);
        integer a =u * u;a.mod2pow(m); 
        w = (1 << (w - 1)) - 1;
        for (int i = 0; i < w; i++)
        {
            integer tmp=a * b.back();
            tmp.mod2pow(m);
            b.push_back(tmp);
        }
        const int ls = l.size();
        a = b[l[ls - 2]];
        for (int j = 0; j < l[ls - 1]; j++) { a = a * a;a.mod2pow(m); }
        for (int i = ls - 4; i > -1; i -= 2)
        {
            a = a * b[l[i]];
            a.mod2pow(m);
            for (int j = 0; j < l[i + 1]; j++)
            {
                a = a * a;
                a.mod2pow(m);
            }
        }
        return a;
    }
public:
    mont a;
    mod(const integer& x)//>1
    {
        if(x.num[0]&1){a1=1,m=0;a.p=x;a.init();return;}
        m=x.ctz();
        integer y=x;y.div2pow(m);
        if(a1=(y.num.size()>1||y.num.back()>1))
        {
            a.p=std::move(y);a.init();
            integer tmp=1;tmp.mul2pow(m);
            key=inv(a.p,tmp);
        }
    }
    pairs in(const integer& x)
    {
        pairs ans;
        if (a1)
        {
            ans.f = a.in(x);
        }
        if (m)
        {
            ans.s = integer::view(x.num.data(),0,std::min((int)x.num.size()-1,m/Blen),x.sign);
            ans.s.mod2pow(m);
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
        if (m)
        {
            ans.s = u.s * v.s;
            ans.s.mod2pow(m);
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
        if (m)
        {
            ans.s = u.s+v.s;
            ans.s.mod2pow(m);
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
        if (m)
        {
            ans.s = u.s-v.s;
            ans.s.mod2pow(m);
        }
        return ans;
    }
    integer out(const pairs& u)//不保证返回正数
    {
        if (!a1) {  return u.s;}
        integer x = a.out(u.f);
        if (!m) { return x; }
        integer tmp=(u.s - x) * key;
        tmp.mod2pow(m);
        return x + a.p * tmp;
    }
    pairs one()
    {
        pairs ans;
        if (a1)
        {
            ans.f = a.f;
        }
        if (m)
        {
            ans.s = 1;
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
        if (m)
        {
            ans.s = 0;
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
        if (m)
        {
            ans.s = pow_exponent_m(x.s, l, w);
        }
        return ans;
    }
};
std::vector<int> exponent(const integer&a, int w = 5)
{
    std::vector<int>l;
    int step = 0, now = 0, cnt = 0;
    for(int i=0;i<a.num.size();i++)
    {
        int m=a.num[i];bool last=(i+1==a.num.size());
        for(int j=0;j<Blen&&(!last||m);j++,m>>=1)
        {
            int x =  m& 1;
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
int getw(const  integer&n)
{
    return std::min((int)n.num.size(), 4) + 1;
}
integer power(const integer& x, const integer& n, const integer& p)//p>1
{
    if (n.num.size()==1)
    {
        int m = n.num.back();
        integer a(x.mod_positive(p)), ans = m & 1 ? a : 1;
        while (m >>= 1)
        {
            a = a * a % p;
            if (m & 1) { ans = ans * a % p; }
        }
        return ans;
    }
    mod q(p); int w = getw(n);
    integer ans=q.out(q.pow_exponent(q.in(x), exponent(n, w), w));
    return ans.sign<0&&ans.num.back()?ans+p:ans;
}
#if more
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
    if (n.num[0] % 2 == 0) { return 0; }
    integer d(n - 1);
    int s = d.ctz();
    d.div2pow(s);
    mont q; q.p = n; q.init();
    const integer flag = q.p - q.f;
    int w = getw(n);
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
    integer x(a.mod_positive(b)), y(b);
    int m = 1;
    while (x.num.back())
    {
        int j = y.num[0] % 8;
        if (x.num[0] % 2 == 0)
        {
            int k = x.ctz();
            x.div2pow(k);
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
    if (x.num.size() == 1 && x.num[0] == 2) { return 1; }
    integer u;
    if (x.num[0] % 2 == 0 ||
        ((inside(x.num[0] % 128, {1, 9, 17, 25, 33, 41, 49, 57, 65, 73, 81, 89, 97, 105, 113, 121}) &&
            ((u = x.fsqrt()) * u).num == x.num))) {
        return 0;
    }
    int t = 5, s = 1;
    while (jacobi(t * s, x) != -1)
    {
        t += 2; s = -s;
    }
    std::vector<int>l = exponent(x + 1, 1);
    mont q; q.p = x; q.init(); u = q.f;
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
    if (x.absbigger(integer(341550071728321), 1)) { return miller(x, { 2,3 }) && lucas(x); }
    return miller(x);
}
integer shanks_sqrt(const integer& a, const integer& p)//p奇质数
{
    int n = 2;
    while (jacobi(n, p) != -1) { n++; }
    integer t = p - 1;
    int e = t.ctz();t.div2pow(e);
    mont q; q.p = p; q.init();
    int w = getw(t);
    integer y = q.pow_exponent(q.in(n), exponent(t,w),w), a0 = q.in(a),
        x = q.pow_exponent(a0, exponent((t - 1) / 2,w),w), b = q.out(a0 * q.out(x * x));
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
    if (n.num.size() > 7) { std::cout << "qs too big"; exit(0); }//qs知乎学的,10^47之内保证不超过10s
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
    std::vector<equa>smooth;
    std::vector<equa>flaw; int fs = 0;
    std::vector<int>index;
    std::vector<int>pcnt;
    int ta = 0, tb = 0, tc = 0; float bound;
    auto insert_lambda =
        [&q, &prime, &smooth, &flaw, &fs, &n, &ok, &ta, &tb, &tc, &bound, &index, &pcnt](const integer& x, integer y) -> integer {
        int  ps = prime.size();
        equa e; e.mask.assign(ps / 64 + (bool)(ps % 64 != 0), 0);
        int z = 0;
        for (; z < ps && (y.num.size() > 1 || y.num[0] > 1); z++)
        {
            integer r, p = prime[z], q = integer::div_native(y, p, r);
            int cnt = 0;
            while (!r.num.back())
            {
                cnt++;
                std::swap(y, q);
                q = integer::div_native(y, p, r);
            }
            pcnt[z] = cnt;
        }
        if (y.num.size() > 1 || (y.num[0]>(Base>>4))) { tb++; return 0; }
        int j = 0;
        if (y.num[0] > 1)
        {
            tb++; for (; j < fs && flaw[j].x.sign != y.num[0]; j++); if (j >= flaw.size() || !miller(y.num[0], { 2 })) { return 0; }
        }
        else { ta++; }
        e.x = q.in(x);
        int  c = 0, r, cend = e.mask.size() * 64;
        for (; z < ps; z++) { pcnt[z] = 0; }
        for (int i = 0; i < index.size(); i++) { pcnt[index[i]]++; }
        e.y = y.sign;
        for (int i = 0; i < pcnt.size(); i++)
        {
            int cnt = pcnt[i], j = cnt & 1, t = cnt >> 1;
            e.mask[i / 64] |= ((unsigned long long)j) << (i % 64);
            if (t) { e.y = e.y * power(prime[i], t) % n; }
        }
        if (y.sign < 0) { e.y = e.y + n; }
        if (y.num[0] > 1)
        {
            equa& s = flaw[j];
            if (s.x.sign < 0) { e.x.sign = y.num[0]; flaw[j] = e; fs++; return 0; }
            tc++;
            s.x.sign = 1; e.x = q.out(e.x * s.x); s.x.sign = y.num[0];
            e.y = (e.y * y.num[0]) * s.y;
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
                    if (u & 1)
                    {
                        l.push_back(prime[c + r]);
                    }
                }
                e.mask[c / 64] ^= si;
            }
            e.y = mul(l)% n * (e.y * s.y) % n;
            for (c = pre; c < cend && !e.mask[c / 64]; c += 64);
            if (c < cend) { for (bit = e.mask[c / 64], r = 0; r < ps - c && r < 64 && !(bit & 1); r++, bit >>= 1); }
        }
        if (index >= ps)
        {
            ok--;
            return gcd(q.out(e.x) + e.y, n);
        }
        smooth[index] = e;
        //std::cout << "\n    relation:" << ta + tc << " need:" << prime.size() << " |from " <<tb << " fail extract " << fs << " flaw and " << tc << " relation\n";
        const int level = 2;
        if (tb > (ta << (level + 1))) { bound += 0.1; }if ((ta << level) > tb) { bound -= 0.1; }//倾向于错误多些,快
        return 0;
        };
    auto getlog = [](const integer& x)->float {
        int l = x.num.size();
        return l == 1 ? log(x.num[0]) : log((ll)x.num.back() * Base + x.num[l - 2]) + (l - 2) * log(Base);
        };
    integer d = n.fsqrt();
    B = pow(getlog(d) - 2, 1.67 + 0.106* n.num.size());
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
    d = d * 14142 / (integer(m) * 10000);
    q.p = n; q.init();
    pcnt.resize(prime.size());
    smooth.assign(prime.size(), equa(n.num.size(), prime.size()));
    flaw.assign(prime.size() << 2, equa(n.num.size(), prime.size()));
    ok = 40;
    bound = getlog(d) + 2 * log(m) - 13;
    integer a, b; int start = 30, gap = 1;
    std::vector<integer>crtbase;
    std::vector<int>inva; inva.resize(root.size());
    auto geta = [&a, &b, &d, &prime, &index, &start, &gap, &root, &crtbase, &inva]()->bool {
        auto inv_fsp = [](const integer& a, int p)->int {
            int q = (a % p).num[0];
            int h = p - 2, inva = h & 1 ? q : 1;
            while (h >>= 1)
            {
                q = q * q % p;
                if (h & 1) { inva = q * inva % p; }
            }
            return inva;
            };
        if (start > std::min((int)prime.size(), 500)) { start = std::min((int)prime.size() / 4, 30); }
        index.clear(); crtbase.clear();
        int i = start;
        integer tmp = prime[i];
        do
        {
            index.push_back(i);
            i += gap + rand() % 8;
            if (i + 1 > prime.size())
            {
                gap++;
                if (gap > prime.size() / 10) { gap = 1; }
                return 0;
            }
            std::swap(a, tmp);
            tmp = a * prime[i];
        } while (d.absbigger(a * prime.back(), 0));
        int j = (d / a).num[0];
        if (prime.size() > 200 && j < 10) { start++; return 0; }//psize太小会卡死
        for (i = 2; i + 1 < prime.size() && prime[i] < j; i++);
        float low = (float)prime[i - 1] / j, high = (float)prime[i] / j;
        if (abs(1 - low) < abs(1 - high)) { i--; }
        for (j = 0; j < index.size() && index[j] != i; j++);
        if (j != index.size()) { start++; return 0; }
        a = a * prime[i];
        index.push_back(i);
        for (i = 0, b = 0; i < index.size(); i++)
        {
            integer p = prime[index[i]], mi = a / p, invi = inv_fsp(mi, p.num[0]), ti = (invi * root[index[i] - 1]) % p * mi % a;
            b = b + ti;
            crtbase.push_back(ti);
        }b = b % a;
        for (int i = 0; i < root.size(); i++)
        {
            inva[i] = inv_fsp(a, prime[i + 1]);
        }
        start++;
        return 1;
        };
    int cnt = 0;
    auto getb = [&a, &b, &prime, &root, &index, &cnt, &crtbase]()->bool {
        cnt++;
        if (cnt + 1 > (1 << (index.size() - 1))) { return 0; }//只要a-b和b中一个
        int u = 1, pos = 0;
        while ((cnt & u) == 0) { u <<= 1; pos++; }
        u = cnt & (u << 1);
        integer& diff = crtbase[pos];
        if (!u) { diff.sign = -1; }
        else { diff.sign = 1; }
        b = (b + (diff + diff)) % a;
        if (b.sign < 0) { b = b + a; }
        return 1;
        };
    while (ok > 0)
    {
        if (!getb()) { cnt = 0; while (!geta()); }
        integer c = (b * b - n) / a;
        std::vector<float>pos(m), neg(m);
        for (int i = 0; i < root.size(); i++)
        {
            int invai = inva[i]; if (!invai) { continue; }
            int p = prime[i + 1], h = (b % p).num[0];
            int x0 = (-h - root[i]) * invai % p;
            int x1 = (-h + root[i]) * invai % p;
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
                g = insert_lambda(b + a * i, (a * ((ll)i * i) + b * (i << 1) + c));
                if (check_lambda()) { return g; }
            }
            if (neg[i] > bound)
            {
                g = insert_lambda(b - a * i, (a * ((ll)i * i) - b * (i << 1) + c));
                if (check_lambda()) { return g; }
            }
        }
    }
    std::cout << "factor"; exit(0);//可能是输入质数幂(qs非法输入)，或者极小概率geta产生重复的a导致，或者ok正常耗尽
}
integer ppow(const integer& x)
{
    integer p;
    if (x.num[0] % 2 == 0){p=2;}
    else
    {
        int w=getw(x);
        std::vector<int>b = exponent(x,w);
        ll a = 2;
        mont q; q.p = x; q.init();
        do
        {
            p = gcd(q.out(q.pow_exponent(q.in(a), b,w)) - a, x);
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
        p = integer::divide(p, l.back(), r);
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
    x= k[0];
    for (int i = 1; i < k.size(); i++)
    {
        l.push_back(k[i] - k[i - 1]);
    }
    mont q; q.p = p; q.init();
    for (int a = 2;; a++)
    {
        integer t = q.in(a), y = q.pow_binary(t, x);
        int s = -1;
        while (y.num != q.f.num)
        {
            if (++s == l.size()) { return a; }
            y = q.out(y * q.pow_binary(t, l[s]));
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
#endif
