
#include "integerb.h"

// ============ 全局表与常量定义 ============

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

std::vector<double> table1 = table(1), table3 = table(3);

std::vector<cd> unit = []() {
	std::vector<cd> l;
	for (int i = 0; i < log2lenmax; i++) {
		double t = -(double)3.1415926535897932384626433832795 / (1 <<(i+1));
		l.push_back(cd(cos(t), sin(t)));
	}
	return l;
}();

cd tb[log2lenmax];
integer::view integer::lastview;
double* integer::lasty;
// ============ c2 的成员函数 ============

void c2::first(const c2& y)
{
	double u = r.a * y.r.a + i.a * y.i.a;
	double v = r.a * y.i.a + i.a * y.r.a;
	r.a = u, i.a = v;
	u = r.b * y.r.b - i.b * y.i.b;
	v = r.b * y.i.b + i.b * y.r.b;
	r.b = u, i.b = v;
}

void c2::second(const c2& y)
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

// ============ fft / ifft / middle ============

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

void middle(c2*xx, c2*yy, int n4)
{
	xx[0].first(yy[0]); xx[1].second(yy[1]);
	int l = 0,m = n4<<1;
	while (!(m & 1)) { m >>= 1, l++; }
	int  lk=l -3;
	int index;
	auto next = [&index, &n4](int& k)->cd {
		int bit = n4, cnt = 0;
		while (k & bit) { k ^= bit, bit >>= 1; cnt++; }k ^= bit;
		index -= cnt; tb[index + 1] = tb[index] * unit[cnt];
		return tb[++index];
	};
	for (int len = 2; len < n4; len <<= 1, lk--)
	{
		int k = 1 << lk; index = 0, tb[0] = unit[l-1 - lk];
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
	}
}
// ============ integer 私有 static 辅助函数 ============
template <bool WantQ>
void integer::quotient(const view& a, const view& b, integer& now)
{
	int m = 1, t = b.ptr[b.len - 1];
	while ((t << m) < Base) { m++; }m -= 1;
	int m_5 = m + 5;
	ull bb = ((ull)b.ptr[b.len - 1] << m_5) + (((ull)b.ptr[b.len - 2] << m) >> (Blen - 5));
	bool needcopy=now.num.data()!=a.ptr||now.num.size()!=a.len;
	now.num.reserve(a.len + 1);
	if(needcopy){now = a;}
	now.num.push_back(0); int ln;
	for (; (t = (ln = now.num.size()) - b.len - 1) > -1; now.num.pop_back())
	{
		int q = (((((ull)now.num.back() << Blen) + now.num[ln - 2]) << m_5) + (((ull)now.num[ln - 3] << m) >> (Blen - 5))) / bb;
		if (!q) { continue; }
		int*start=now.num.data()+t;
		ll k = 0;
		for (int i = 0; i < b.len; i++)
		{
			k +=(ll)start[i]-(ll)b.ptr[i] * q;
			start[i] = k & Bmask, k >>= Blen;
		}
		start[b.len] += k;
		int top=b.len;
		while(top>-1&&start[top]==0){top--;}
		if(top>-1&&start[top]<0)
		{
			int carry=0;
			for(int i=0;i<b.len;i++)
			{
				carry+=start[i]+b.ptr[i];
				start[i]=carry&Bmask;
				carry>>=Blen;
			}q--;
		}
		if constexpr (WantQ)
		{
			num[t] = q;
		}
	}
	while (now.num.back() == 0 && now.num.size() > 1) { now.num.pop_back(); }
}


void integer::div_short(int*anum,int alen,const int*bnum,int blen,integer&q,int qstart)
{
	int l=blen+1;
	ull bb = ((ull)bnum[blen - 1] <<5) + (((ull)bnum[blen - 2]) >> (Blen - 5));
	for(int qindex=alen-blen-1;qindex>-1;qindex--)
	{
		int q1 = (((((ull)anum[blen+qindex]<< Blen) +anum[blen+qindex-1]) <<5) + (((ull)anum[blen+qindex-2]) >> (Blen - 5))) / bb;
		if(!q1){continue;}
		int*start=anum+qindex;
		ll k=0;
		for(int i=0;i<blen;i++)
		{
			k+=(ll)start[i]-(ll)q1*bnum[i];
			start[i]=k&Bmask;
			k>>=Blen;
		}
		start[blen]+=k;
		int top=blen;
		while(top>-1&&start[top]==0){top--;}
		if(top>-1&&start[top]<0)
		{
			int carry=0;
			for(int i=0;i<blen;i++)
			{
				carry+=start[i]+bnum[i];
				start[i]=carry&Bmask;
				carry>>=Blen;
			}q1--;
			start[blen]=0;
		}
		q.num[qstart+qindex]=q1;
	}		
}

void integer::div_4n_2n(int*anum,int alen,const int*bnum,int blen,integer&q,int qstart)
{
	int qlen=alen-blen;
	if(qlen<15||blen<15){div_short(anum,alen,bnum,blen,q,qstart);return;}
	int qlen1=qlen>>1,alen1=blen+qlen1,ashift1=alen-alen1;
	div_3n_2n(anum+ashift1,alen1,bnum,blen,q,qstart+ashift1);
	alen-=qlen1;
	int j=alen;
	while(j>1&&anum[j-1]==0){j--;}
	if(j!=alen){alen=j+(anum[j-1]>=bnum[blen-1]);}
	div_3n_2n(anum,alen,bnum,blen,q,qstart);
}

bool integer::qlen_overflow(int*anum,int alen,const int*bnum,int blen,integer&q,int qstart)
{
	int i=blen-1,qlen=alen-blen;
	while(i>-1&&anum[i+qlen]==bnum[i]){i--;}
	if(i<0)
	{
		for(i=0;i<qlen;i++)
		{
			q.num[qstart+i]=Bmask;
		}
		return 1;
	}
	return 0;
}

void integer::div_3n_2n(int*anum,int alen,const int*bnum,int blen,integer&q,int qstart)
{
	int qlen=alen-blen;
	if(qlen<15||blen<15){div_short(anum,alen,bnum,blen,q,qstart);return;}
	int shift=blen-qlen;
	bool flag=qlen_overflow(anum+shift,alen-shift,bnum+shift,blen-shift,q,qstart);
	if(!flag){div_4n_2n(anum+shift,alen-shift,bnum+shift,blen-shift,q,qstart);}
	view qview(q.num.data(),qstart,qstart+qlen-1,1);
	integer prod=karamul(qview,view(bnum,0,(flag?blen:shift)-1,1));
	int cnt=0,s=abssub(anum,alen,prod.num.data(),prod.num.size());
	while(s<0)
	{
		s=-abssub(anum,alen,bnum,blen);
		if(cnt++>2){std::cout<<"div_3n_2n ";exit(0);}
	}
	if(cnt){abssub(q.num.data()+qstart,qlen,&cnt,1);}
}
// ============ Fermat / NTT 相关私有函数 ============

 void integer::transform2_fermat(int*a,int*b,int m)
    {
        bool bm=b[m];
        if(a[m]&&bm)
        {
            for(int i=0;i<m;i++){a[i]=Bmask;}
            a[m]=b[m]=0;
            return;
        }
        int kadd=0;
        bool ksub=0;
        for(int i=0;i<m;i++)
        {
            kadd+=a[i]+b[i];
            b[i]=a[i]-b[i]-ksub;
            a[i]=kadd&Bmask;
            kadd>>=Blen;
            ksub=b[i]<0;
            b[i]+=ksub<<Blen;
        }
        b[m]=a[m]-b[m]-ksub;
        if(kadd||bm||a[m])
        {
            int i=0;
            for(;i<m&&a[i]==0;i++);
            if(i!=m)
            {
                for(int j=0;j<i;j++){a[j]=Bmask;}
                a[i]--,a[m]=0;
            }
            else{a[m]=1;}
        }
        if(b[m]<0)
        {
            b[m]=0;
            int j=0;
            for(;j<m&&b[j]==Bmask;j++){b[j]=0;}
            b[j]++;
        }    
    }
    void  integer::neg_fermat(int*a,int*out,int m)
    {
        if(a[m])
        {
            out[0]=1;
            for(int i=1;i<=m;i++){out[i]=0;}
            return;
        }
        int k=2;
        for(int i=0;i<m;i++)
        {
            k+=Bmask-a[i];
            out[i]=k&Bmask;
            k>>=Blen;
        }
        out[m]=k;
        if(out[m]&&out[0])
        {
            out[m]=out[0]=0;
        }
    }
    void integer::mulBasepowmod_fermat(int*a,int q,int r,int m,int*tmp)
    {
        if(!r){if(q&1){neg_fermat(a,a,m);}return;}
        if(q&1)
        {
            neg_fermat(a,tmp,m);
        }
        else
        {
            for(int i=0;i<=m;i++){tmp[i]=a[i];}
        }
        int i=0;
        for(;i<r&&tmp[i+m-r]==0;i++){a[i]=0;}
        bool k=(i!=r);
        if(k)
        {
            a[i]=Base-tmp[i+m-r];i++;
            for(;i<r;i++){a[i]=Bmask-tmp[i+m-r];}
        }
        a[r]=tmp[0]-tmp[m]-k;
        k=a[r]<0;
        a[r]+=k<<Blen;i++;
        if(k)
        {
            for(;i<m&&tmp[i-r]==0;i++){a[i]=Bmask;}
            if(i==m)
            {
                a[m]=0;
                int j=0;
                for(;j<m&&a[j]==Bmask;j++){a[j]=0;}
                a[j]++;
                return;
            }
            a[i]=tmp[i-r]-1;i++;
        }
        for(;i<m;i++){a[i]=tmp[i-r];}
        a[i]=0;
    }
    void integer::mod_fermat(int*in,int*out,int m,int lb)
    {
        int*a=in,*b=in+m;
        bool k=0;
        int i=0;
        for(;i<lb;i++)
        {
            out[i]=a[i]-b[i]-k;
            k=out[i]<0;
            out[i]+=k<<Blen;
        }
        if(k)
        {
            for(;i<m&&a[i]==0;i++){out[i]=Bmask;}
            if(i==m)
            {
                out[m]=0;
                int j=0;
                for(;j<m&&out[j]==Bmask;j++){out[j]=0;}
                out[j]++;
                return;
            }
            out[i]=a[i]-1;i++;
        }
        for(;i<m;i++){out[i]=a[i];}
        out[i]=0;
    }   
    void integer::ntt_dfs(int*x,int mlog2,int lenlog2,int*tmp)
    {
        int m=1<<mlog2,l=m+1;
        if(lenlog2==1){transform2_fermat(x,x+l,m);return;}
        int half=l<<(lenlog2-1);
        for(int k=0,end=1<<(lenlog2-1),curexp=0,expstep=1<<(1+mlog2-lenlog2);k<end;k++)
        {
            transform2_fermat(x,x+half,m);
            mulBasepowmod_fermat(x+half,curexp>>mlog2,curexp&(m-1),m,tmp);
            x+=l,curexp+=expstep;
        }
        ntt_dfs(x-half,mlog2,lenlog2-1,tmp);
        ntt_dfs(x,mlog2,lenlog2-1,tmp);
    }
    void  integer::intt_dfs(int*x,int mlog2,int lenlog2,int*tmp)
    {
        int m=1<<mlog2,l=m+1;
        if(lenlog2==1){transform2_fermat(x,x+l,m);return;}
        int half=l<<(lenlog2-1);
        intt_dfs(x,mlog2,lenlog2-1,tmp);
        intt_dfs(x+half,mlog2,lenlog2-1,tmp);
        for(int k=0,end=1<<(lenlog2-1),curexp=m<<1,expstep=1<<(1+mlog2-lenlog2);k<end;k++)
        {
            mulBasepowmod_fermat(x+half,curexp>>mlog2,curexp&(m-1),m,tmp);
            transform2_fermat(x,x+half,m);
            x+=l,curexp-=expstep;
        }
    }
    void integer::invlen_fermat(int*x,int mlog2,int lenlog2,int*tmp)
    {
        int m=1<<mlog2,l=m+1,n0=Blen-lenlog2,q1=2*m-1,r1=q1&(m-1);q1>>=mlog2;
        for(int i=0,len=1<<lenlog2;i<len;i++,x+=l)
        {
            int m0=Blen-n0,mask=(1<<m0)-1;
            tmp[m+1]=x[m]>>m0;
            for(int j=m;j>0;j--)
            {
                int high=mask&x[j],low=x[j-1]>>m0;
                tmp[j]=low|(high<<n0);
            }
            tmp[0]=(x[0]&mask)<<n0;
            mod_fermat(tmp,x,m,2);
            mulBasepowmod_fermat(x,q1,r1,m,tmp);
        }
    }
    int* integer::mk(const view&a,int k,int m,int len)
    {
        int l=m+1;
        int*x=new int[len*(m+1)]();
        int i=0,j=0;
        for(;i+k<=a.len;i+=k,j+=l)
        {
            for(int t=0;t<k;t++)
            {
                x[j+t]=a.ptr[i+t];
            }
        }
        if(i!=a.len)
        {
            int k0=a.len-i;
            for(int t=0;t<k0;t++)
            {
                x[j+t]=a.ptr[i+t];
            }
        }
        return x;
    }
    cd integer::Wi(int i,int&index,int&ilog2,int mlog2)//w(i,4m),i++从0到m-1
    {
        if(!i){return 1;}
        if(!(i&(i-1)))
        {
            index=0;
            tb[0]=unit[mlog2-ilog2];
            ilog2++;
            return tb[0];
        }
        i--;
        int bit=1,cnt=0;
        while(i&bit){i^=bit,bit<<=1,cnt++;}i^=bit;
        index-=cnt;
        tb[index+1]=tb[index]*unit[mlog2-cnt];
        return tb[++index];
    }
    void integer::b14rrii(int*l,int mlog2,c2*xx)
    {
        int m2=1<<(mlog2-1),i=0,index,ilog2=0;
        for(int i=0;i<m2;i++)
        {
            cd u=conj(Wi(2*i,index,ilog2,mlog2))*cd(l[i]&Fmask,l[i+m2]&Fmask);
            cd v=conj(Wi(2*i+1,index,ilog2,mlog2))*cd(l[i]>>Flen,l[i+m2]>>Flen);
            xx[i]=c2(u,v);
        }
    }


double* integer::b14rrii(const view& a, int n)
{
	double* y = new double[n]();
	int i = 0;
	for (; i + 1 < a.len; i += 2)
	{
		y[2 * i] = a.ptr[i] & Fmask;
		y[2 * i + 2] = a.ptr[i] >> Flen;
		y[2 * i + 1] = a.ptr[i + 1] & Fmask;
		y[2 * i + 3] = a.ptr[i + 1] >> Flen;
	}
	if (i != a.len)
	{
		y[2 * i] = a.ptr[i] & Fmask;
		y[2 * i + 2] = a.ptr[i] >> Flen;
	}
	return y;
}
// ============ integer 构造 / 赋值 / 字符串转换 / 基础工具 ============

integer::integer(integer&& other) noexcept : sign(other.sign), num(std::move(other.num)) {}

integer::integer(const integer& other) : sign(other.sign), num(other.num) {}

integer& integer::operator=(integer&& other) noexcept
{
	sign = other.sign;
	num = std::move(other.num);
	return *this;
}

integer& integer::operator=(const integer& other)
{
	sign = other.sign;
	num = other.num;
	return *this;
}

integer::integer(const view& v)
{
	num.assign(v.ptr, v.ptr + v.len);
	sign = v.sign;
}

integer::integer(ll n)
{
	if (n < 0) { n = -n; sign = -1; }
	do
	{
		num.push_back(n & Bmask);
		n >>= Blen;
	} while (n);
}

void integer::get()
{
	std::string s;
	std::cin >> s;
	*this = integer(s);
}

integer::integer(const std::string& s)
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
	e8tob28();
}

void integer::e8tob28()
{
	while(num.size()>1&&num.back()==0){num.pop_back();}
	int n = 1; while (n < num.size()) { n <<= 1; }
	num.resize(n);
	integer h; h.num.push_back(e8);
	for (int len = 1; len < n;)
	{
		int len2 = len << 1;
		for (int i = 0; i < n; i += len2)
		{
			view high(num.data() + i, len, len2 - 1, 1);
			view low(num.data() + i, 0, len - 1, 1);
			integer tmp = karamul(high, h);
			tmp = addorsub(tmp.num.data(), tmp.num.size(), 1, low.ptr, low.len, 1, 1);
			int j = 0;
			for (; j < tmp.num.size(); j++)
			{
				num[i + j] = tmp.num[j];
			}
			for (; j < len2; j++) { num[i + j] = 0; }
		}
		len = len2;
		if (len < n) { h = h * h; }
		else { break; }
	}
	while (num.back() == 0 && num.size() > 1) { num.pop_back(); }
}

ll integer::toll()const
{
	int ln = num.size();
	if (ln > 2) { std::cout << "toll"; exit(0); }
	ll a = num.back();
	if (ln > 1) { a = (a << Blen) + num[0]; }
	return sign * a;
}

double integer::getlog()const
{
	int l = num.size();
	return l == 1 ? log(num[0]) : log((ll)num.back() * Base + num[l - 2]) + (l - 2) * log(Base);
};

std::string integer::tostring(int m,int floatbe8)const
{
	if (m < 1) { m = num.size() << 2; }
	std::string s;
	if (sign != 1) { s.push_back('-'); }
	if (num.back()) {
		integer E8; E8.num.push_back(e8);
		int l = num.size(), n = 0; integer tmp;
		if (l > m)
		{
			n = (l - m) * (log(Base) / log(e8)) - 1;
		}
		if (n < 1)
		{
			n = 0;
			tmp = *this;
		}
		else { tmp = *this / shiftpow(E8, n, num.size(), l); }
		if (tmp.absbigger(E8, 1))
		{
			std::vector<integer>p; p.push_back(E8);
			int t = 2;
			while (tmp.absbigger(p.back(), 0) && p.back().num.size() * 2 - 1 <= tmp.num.size())
			{
				t <<= 1; p.push_back(p.back() * p.back());
			}
			if (p.back().absbigger(tmp, 0)) { p.pop_back(); t /= 2; }
			tmp.num.resize(t);
			for (int len = t; len > 1;)
			{
				int len2 = len >> 1;
				for (int i = 0; i < t; i += len2)
				{
					view a(tmp.num.data() + i, 0, len - 1, 1);
					integer r, q = divide(a, p.back(), r);
					int j = 0;
					for (; j < r.num.size(); j++)
					{
						tmp.num[i + j] = r.num[j];
					}
					for (; j < len2; j++) { tmp.num[i + j] = 0; }
					j = 0, i += len2;
					for (; j < q.num.size(); j++)
					{
						tmp.num[i + j] = q.num[j];
					}
					for (; j < len2; j++) { tmp.num[i + j] = 0; }
				}
				len = len2;
				if (len > 1) { p.pop_back(); }
				else { break; }
			}
			while (tmp.num.back() == 0) { tmp.num.pop_back(); }
		}
		s.reserve(tmp.num.size() * 8 + 20);
		int cnt = 0; int a[8];
		a[7] = tmp.num.back();
		for (int j = 7; j > 0; j--)
		{
			a[j - 1] = a[j] / 10;
		}
		for (; cnt < 7 && a[cnt] == 0; cnt++);
		if (n||floatbe8) {
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
		if (n||floatbe8)
		{
			s.push_back('e');
			s += std::to_string(((int)tmp.num.size() + n+floatbe8) * 8 - cnt);
		}
	}
	else {
		s.push_back('0');
	}
	return s;
}

void integer::print(bool stdendl, int m)const
{
	std::cout << tostring(m);
	if (stdendl) { std::cout << '\n'; }
	else { std::cout << ' '; }
}

bool integer::absbigger(const view& that, bool equal) const
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

integer integer::operator-()const
{
	integer a = *this; a.sign = -sign;
	return a;
}

void integer::shiftadd(const view& b, int i)
{
	if(b.len+i>num.size()){std::cout<<"shiftadd";exit(0);}
	int lb = b.len, k = 0;
	for (int j = 0; j < lb; j++)
	{
		num[j + i] += b.ptr[j] + k;
		if (num[i + j] >= Base) { num[i + j] -= Base; k = 1; }
		else { k = 0; }
	}
	if (k) { for (int j = i + lb; j < num.size() && ++num[j] == Base; j++) { num[j] = 0; } }
}
// ============ integer 加减乘除运算符 + abssub / addorsub / div_native / multiply ============

integer integer::operator+(const integer& that)const
{
	return addorsub(num.data(), num.size(), sign, that.num.data(), that.num.size(), that.sign, 1);
}

integer integer::operator-(const integer& that)const
{
	return addorsub(num.data(), num.size(), sign, that.num.data(), that.num.size(), that.sign, 0);
}

integer integer::operator*(const integer& that)const
{
	return karamul(*this, that);
}

integer integer::operator/(const integer& that)const
{
	integer r;
	return divide(*this, that, r);
}

integer integer::operator%(const integer& that)const
{
	integer r;
	divide(*this, that, r,0);
	return r;
}

void integer::operator%=(const integer& that)
{
	divide(*this, that, *this,0);
}
template <bool WantQ>
integer integer::div_native(const view& a, const view& b, integer& r)
{
	if (!b.ptr[b.len - 1]) { std::cout << "div0"; exit(0); }
	if(a.len<b.len){r=a;return 0;}
	integer c;
	if constexpr (WantQ)
	{
		c.num.resize(std::max(a.len - b.len + 1, 1));
		c.sign = a.sign * b.sign;
	}
	if (b.len == 1)
	{
		ull n = b.ptr[0], now = 0;
		for (int i = a.len - 1; i > -1; i--)
		{
			now = (now << Blen) + a.ptr[i];
			if constexpr (WantQ)
			{
				c.num[i] = now / n;
			}
			now %= n;
		}
		r.num.resize(1),r.num[0]= now;
	}
	else
	{
		c.quotient<WantQ>(a, b, r);
	}
	while (c.num.size() > 1 && c.num.back() == 0) { c.num.pop_back(); }
	r.sign = a.sign;
	return c;
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
		result.num.push_back(0);
		result.shiftadd(view(b,0,lb-1,1),0);
		if(!result.num.back()){result.num.pop_back();}
	}
	else
	{
		if(abssub(result.num.data(), result.num.size(), b, lb)<0){result.sign*=-1;}
		while (result.num.back() == 0 && result.num.size() > 1) { result.num.pop_back(); }
	}
	if (!add) { result.sign *= change; }
	return result;
}

int integer::abssub(int* a, int la, const int* b, int lb,int bshift)//a-(b<<bshift)
{
	if (la < lb+bshift) { std::cout << "abssub"; exit(0); }
	for (int i = 0; i < lb; i++)
	{
		a[i+bshift] -= b[i];
	}
	while (a[la - 1] == 0 && la > 1) { la--; }
	if(a[la-1]==0){return 0;}
	int s; bool k = 0;
	if (a[la - 1] < 0)
	{
		s = -1;
		for (int i = 0; i < la; i++)
		{
			a[i] = -a[i] - k;
			k = a[i] < 0;
			a[i] += k << Blen;
		}
	}
	else
	{
		int i = bshift; s = 1;
		for (int end=std::min(la,lb+bshift); i <end; i++)
		{
			a[i] -= k;
			k = a[i] < 0;
			a[i] += k << Blen;
		}
		for (; i < la && k; i++)
		{
			a[i] -= k;
			k = a[i] < 0;
			a[i] += k << Blen;
		}
	}
	return s;
}

integer integer::multiply(const int* a, int la, const int* b, int lb, int sign,int i,int iend)
{
	if (!lb || (la + a[0] < 2) || (lb + b[0] < 2)) { return 0; }
	int l = la+lb-1;
	integer c;
	c.sign = sign;
	c.num.resize(l + 1);
	if(l>3&&a==b&&la==lb){square(a,la,c,i,iend);return c;}
	ull k = 0;
	int n = la - 1, m = lb - 1;
	for (int end=std::min(m,iend); i < end; i++)
	{
		int j = i;
		for (; j > 3; j -= 4)
		{
			k += dot4(a+i-j,b+j);
		}
		for (; j > -1; j--)
		{
			k += (ll)a[i - j] * b[j];
		}
		c.num[i]+= k & Bmask;
		k >>= Blen;
	}
	for (int end=std::min(n,iend); i < end; i++)
	{
		int j = m;
		for (; j > 3; j -= 4)
		{
			k += dot4(a+i-j,b+j);
		}
		for (; j > -1; j--)
		{
			k += (ll)a[i - j] * b[j];
		}
		c.num[i]+= k & Bmask;
		k >>= Blen;
	}
	for (int end=std::min(l,iend); i < end; i++)
	{
		int j = m;
		for (; j > i - n + 3; j -= 4)
		{
			k += dot4(a+i-j,b+j);
		}
		for (; j > i - n - 1; j--)
		{
			k += (ll)a[i - j] * b[j];
		}
		c.num[i]+= k & Bmask;
		k >>= Blen;
	}
	if (k) { c.num[l] = k; }
	else { c.num.pop_back(); }
	return c;
}
void integer::square(const int* a, int la, integer& c, int i, int iend)
{
    ull k = 0;int icopy=i;
    for(int end = std::min(la-1, iend); i < end; i++)
    {
        int j = 0, i2 = (i>>1)+(i&1);
        for(; j + 3 < i2; j += 4)
        {
            k+= dot4(a+j,a+i-j)<<1;
        }
        for(; j < i2; j++)
        {
            k += ((ll)a[i-j] * a[j]) << 1;
        }
        c.num[i] = k & Bmask;
        k >>= Blen;
    }
    int l = (la<<1) - 1;
    for(int end = std::min(l, iend); i < end; i++)
    {
        int j = i - la + 1, i2  = (i>>1)+(i&1);
        for(; j + 3 < i2; j += 4)
        {
            k+= dot4(a+j,a+i-j)<<1;
        }
        for(; j < i2; j++)
        {
            k += ((ll)a[i-j] * a[j]) << 1;
        }
        c.num[i] = k & Bmask;
        k >>= Blen;
    }
    c.num[l] = k;
	k=0,i=icopy-(icopy&1);
    for(int j=i>>1,end=std::min(l,iend);i<end;j++,i+=2)
	{
		k+=c.num[i]+(ll)a[j]*a[j];
		c.num[i]=k&Bmask;
		k=c.num[i+1]+(k>>Blen);
		c.num[i+1]=k&Bmask;
		k>>=Blen;
	}
	if(!c.num.back()){c.num.pop_back();}
}
// ============ 乘法核心：karamul / shiftmul / ceil2pow / fftmul ============
integer integer::karamul(view  a, view b,int i,int iend)
{
    bool change=a.len < b.len;
	if (change) { std::swap(a, b); }
	bool canfft=a.len + b.len <= halflen;
	if (b.len < 80 || (b.len < 160 &&((!canfft)||a.len * b.len < 30000))) { return multiply(a.ptr, a.len, b.ptr, b.len, a.sign * b.sign,i,iend); }
	if (canfft) {if (change) { std::swap(a, b); }return fftmul(a, b); }
	if(b.len>halflen*0.7){return nttmul(a,b);}
	int n = (a.len + 1) / 2;
	integer ans;
	if (b.len <= n*1.3)
	{
		int la=halflen-1-b.len;
		int k=a.len/la+(bool)(a.len%la);
		la=a.len/k+bool(a.len%k);
        view lv=lastview;double*ly=lasty;
		trycatch();
		ans.num.resize(a.len+b.len);
		for(int i=0;i<a.len;i+=la)
		{
			view anow(a.ptr,i,std::min(i+la,a.len)-1,1);
			ans.shiftadd(fftmul(anow,b),i);
		}	
		endcatch();
        lastview=lv,lasty=ly;
	}
	else
	{
        int state=lastview.sign;
        lastview.sign=-1;
		view b1(b.ptr, n, b.len - 1, 1);
		view b0(b.ptr, 0, n - 1, 1);
		view a1(a.ptr, n, a.len - 1, 1);
		view a0(a.ptr, 0, n - 1, 1);
		integer c1 = karamul(a1, b1);
		ans=karamul(a0, b0);
		integer tmp=addorsub(a0.ptr, a0.len, 1, a1.ptr, a1.len, 1, 1);
		tmp = karamul(tmp,a.ptr==b.ptr&&a.len==b.len?tmp:addorsub(b0.ptr, b0.len, 1, b1.ptr, b1.len, 1, 1));
		abssub(tmp.num.data(), tmp.num.size(), ans.num.data(),ans.num.size());
		abssub(tmp.num.data(), tmp.num.size(), c1.num.data(), c1.num.size());
		while (tmp.num.back() == 0 && tmp.num.size() > 1) { tmp.num.pop_back(); }
		ans.num.resize(a.len + b.len);
		ans.shiftadd(tmp, n);
		ans.shiftadd(c1, 2 * n);
        lastview.sign=state;
	}
	while (ans.num.size() > 1 && ans.num.back() == 0) ans.num.pop_back();
	ans.sign = a.sign * b.sign;
	return ans;
}

integer integer::shiftmul(const view& a, const view& b, int need, int& offset)
{
	if (need < 1) { std::cout << "shiftmul"; exit(0); }
	int b1 = std::max(0, a.len - need);
	int b2 = std::max(0, b.len - need);
	offset += b1 + b2;
	view aview(a, b1), bview(b, b2);
	return karamul(aview, bview,std::max(0,aview.len+bview.len-need-3));
}
void integer::trycatch()
{
    lastview.sign=0;
}
void integer::endcatch()
{
    if(lastview.sign>0){delete[]lasty;}
    lastview.sign=-1;
}
integer integer::fftmul(const view&a,const view&b)
{
    int m = (a.len+b.len) << 1, n = 256;
	while (n < m) { n <<= 1; }
	if (n > lenmax || n < 0) { std::cout << "fftmul"; exit(0); }
	int n4 = n >> 2;
	double* x0 = b14rrii(a, n), * y0;
	c2* xx = reinterpret_cast<c2*>(&x0[0]), * yy;
	fft(xx, n4);
	const bool same = a.ptr == b.ptr && a.len == b.len;
    const bool uselast=lastview.sign==n&&b.ptr==lastview.ptr&&b.len==lastview.len;
	if (same) { y0 = x0, yy = xx; }
    else if(uselast){y0 = lasty, yy = reinterpret_cast<c2*>(&y0[0]);}
	else { y0 = b14rrii(b, n), yy = reinterpret_cast<c2*>(&y0[0]); fft(yy, n4);}
	middle(xx,yy,n4);
    if(lastview.sign==0){lastview=b;lastview.sign=n;lasty=y0;}
	else if (!same&&!uselast) { delete[]y0; }
	ifft(xx, n4);
	integer c;
	c.num.reserve(m=a.len + b.len + 1); m>>=1;
	c.sign = a.sign * b.sign;
	ll k = 0;
	for (int i = 0; i < m; i++)
	{
		ll t = (ll)(xx[i].r.a + 0.5) + k;
		ll pre = t & Fmask;
		k = t >> Flen;
		t = (ll)(xx[i].i.a + 0.5) + k;
		k = t >> Flen;
		c.num.push_back(((t & Fmask) << Flen) + pre);
		t = (ll)(xx[i].r.b + 0.5) + k;
		pre = t & Fmask;
		k = t >> Flen;
		t = (ll)(xx[i].i.b + 0.5) + k;
		k = t >> Flen;
		c.num.push_back(((t & Fmask) << Flen) + pre);
	}
	delete[]x0;
	while (c.num.size() > 1 && c.num.back() == 0) c.num.pop_back();
	return c;
}
// ============ nttmul（独占一批） ============

integer integer::nttmul(const view&a,const view&b)
    {
        int l=a.len+b.len;
        if(l>>26){std::cout<<"nttmul reject";exit(0);}
        int mlog2=std::min(log2lenmax-1,15),m=1<<mlog2,lenlog2=1;
        while(m>=l&&m>8){m>>=1;mlog2--;}
        int k=(m>>1)-1;
        while(l>(k<<lenlog2))
        {
            lenlog2++;
        }
        int len=1<<lenlog2;
        if(mlog2+1<lenlog2){std::cout<<"nttmul need bigger log2lenmax";exit(0);}
        const bool same=a.ptr==b.ptr&&a.len==b.len;
        int*x=mk(a,k,m,len),*y;
		int*tmp=new int[m+2];
        ntt_dfs(x,mlog2,lenlog2,tmp);
        if(same){y=x;}
        else{y=mk(b,k,m,len);ntt_dfs(y,mlog2,lenlog2,tmp);}
        int n=m<<1,n4=n>>2;
        double*xx=new double[n],*yy;
        c2*xxx=reinterpret_cast<c2*>(&xx[0]),*yyy;
        if(same){yy=xx,yyy=xxx;}
        else{yy=new double[n];yyy=reinterpret_cast<c2*>(&yy[0]);}
        for(int i=0,gap=m+1,end=gap<<lenlog2;i<end;i+=gap)
        {
            if(y[i+m]){neg_fermat(x+i,x+i,m);continue;}
            if(x[i+m]){neg_fermat(y+i,x+i,m);continue;}
            b14rrii(x+i,mlog2,xxx);
            fft(xxx,n4);
            if(!same){b14rrii(y+i,mlog2,yyy);fft(yyy,n4);}
            for(int i=0;i<n4;i++){xxx[i]=xxx[i].mul(yyy[i]);}
            ifft(xxx,n4);
            int index,j=0,jlog2=0;ull k=0;
            for(;j<n4;j++)
            {
                cd w=Wi(2*j,index,jlog2,mlog2);
                xxx[j]=xxx[j].mul(c2(w,Wi(2*j+1,index,jlog2,mlog2)));
                ll t=std::llround(xxx[j].r.a)+k;
                ll pre=t&Fmask;
                k=t>>Flen;
                t=std::llround(xxx[j].r.b)+k;
                k=t>>Flen;
                x[i+j]=((t&Fmask)<<Flen)+pre;
            }
            for(j=0;j<n4;j++)
            {
                ll t=std::llround(xxx[j].i.a)+k;
                ll pre=t&Fmask;
                k=t>>Flen;
                t=std::llround(xxx[j].i.b)+k;
                k=t>>Flen;
                x[i+n4+j]=((t&Fmask)<<Flen)+pre;
            }
            x[i+m]=0;
            int carry=-k;
            for(j=0;j<m&&carry;j++)
            {
                carry+=x[i+j];
                x[i+j]=carry&Bmask;
                carry>>=Blen;
            }
            if(j==m)
            {
                if(carry>0)
                {
                    for(j=0;j<m&&x[i+j]==0;j++);
                    if(j!=m)
                    {
                        for(int t=0;t<j;t++){x[i+t]=Bmask;}
                        x[i+j]--;
                    }
                    else{x[i+m]=1;}
                }
                else
                {
                    for(j=0;j<m&&x[i+j]==Bmask;j++);
                    for(int t=0;t<j;t++){x[i+t]=0;}
                    x[i+j]++;
                }
            }
        }
        delete[]xx;
        if(!same){delete[]yy;delete[]y;}
        intt_dfs(x,mlog2,lenlog2,tmp);
        invlen_fermat(x,mlog2,lenlog2,tmp);
        delete[]tmp;
        integer c;c.sign=a.sign*b.sign;
        c.num.resize(l);
        for(int i=0,gap=m+1,gap1=k;i<len;i++)
        {
            int convi=gap1*i,xi=gap*i;
            for(int j=0;j<gap&&j<l-convi;j++)
            {
                c.num[convi+j]+=x[xi+j];
            }
        }
        delete[]x;
        int carry=0;
        for(int i=0;i<l;i++)
        {
            carry+=c.num[i];
            c.num[i]=carry&Bmask;
            carry>>=Blen;
        }
        while(!c.num.back()&&c.num.size()>1){c.num.pop_back();}
        return c;
    }

// ============ reciprocal / isqrt / root / shiftpow / div_newton / div_bz / divide ============

integer integer::reciprocal(const view& n,int l)
{
	l++;
	int k = 1;
	int t = n.len + k;
	int b;
	integer tmp,xt = div_native<1>(integer(1).shift(3), view(n.ptr, n.len - 2, n.len - 1, n.sign), tmp);
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
		view tmpview(tmp, -(b + k - t));
		abssub(xt.num.data(), xt.num.size(), tmpview.ptr, tmpview.len);
		while (xt.num.back() == 0 && xt.num.size() > 1) { xt.num.pop_back(); }
		t += k;
		k = xt.num.size() - 1;
		while (!xtsize.empty() && xt.num.size() >= xtsize.back()) { xtsize.pop_back(); }
	}
	xt.div2pow(Blen);
	return xt;
}

integer integer::isqrt(bool floor)const
{
	int ns = num.size();
	if (ns < 3)
	{
		ll a = toll(), b = sqrt(a);
		if (ns < 2) { return b; }
		b = (b + 1 + a / (b + 1)) / 2;
		while (b * b > a) { b--; }
		return b;
	}
	int k = 2;
	int f = ns / 2 + (ns % 2 != 0);
	const bool useyk = f > 100;
	int need, b;
	integer xt;
	double guess = (double)((ll)num.back() * Base + num[ns - 2]) * Base + num[ns- 3];
	if (ns % 2)
	{
		xt = sqrt(guess) + 1;
	}
	else
	{
		xt = sqrt(guess * Base + num[ns - 4]) + 1000;
		integer q; need = k + 1;
		do
		{
			b = 0;
			integer tmp = shiftmul(xt, xt, need, b).shift(b);
			view nview(*this, 2 * (f - k));
			abssub(tmp.num.data(), tmp.num.size(), nview.ptr, nview.len);
			while (tmp.num.size() > 1 && tmp.num.back() == 0) { tmp.num.pop_back(); }
			q = tmp / xt; q.div2pow(1);
		} while (q.sign > 0 && q.num.back() && (xt = xt - q).num.back());
	}
	int p = 1;
	integer yt; int l, yk = 0;
	std::vector<int>kplan;
	int k1 = f; k1 += ((k1 & 1) == 0);
	while (k1 > 5) { k1 = k1 / 2 + 1; k1 += ((k1 & 1) == 0); kplan.push_back(k1); }
	while (k < f)
	{
		if (p + k > f) { p = f - k; }
		else if (!kplan.empty() && kplan.back() > k && kplan.back() - k < p)
		{
			p = kplan.back() - k; kplan.pop_back();
		}
		need = xt.num.size() + p + 2, b = 0;
		integer tmp = shiftmul(xt, xt, need, b).shift(p + b);
		view nview(*this, 2 * (f - k) - p);
		int check = abssub(tmp.num.data(), tmp.num.size(), nview.ptr, nview.len);
		if (check < 0) { std::cout << "isqrt"; exit(0); }
		while (tmp.num.size() > 1 && tmp.num.back() == 0) { tmp.num.pop_back(); }
		integer q; bool flag = useyk;
		if (useyk)
		{
			if (yt.num.empty())
			{
				l = 2 * xt.num.size();
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
			if (need > 7)
			{
				q = shiftmul(tmp, yt, need, b);
				view qview(q, l - b);
				integer prod = karamul(qview, xt);
				tmp.sign = abssub(tmp.num.data(), tmp.num.size(), prod.num.data(), prod.num.size());
				while (tmp.num.back() == 0 && tmp.num.size() > 1) { tmp.num.pop_back(); }
				integer ql= divide(tmp, xt, tmp);
				q = addorsub(qview.ptr, qview.len, qview.sign, ql.num.data(), ql.num.size(), ql.sign, 1);
				if (tmp.sign == -1 && tmp.num.back()) {q.addsmall(-1); }
			}
			else { flag = 0; }
		}
		if (!flag)
		{
			q = divide(tmp,xt,tmp);
		}
		q.div2pow(1);
		xt = xt.shift(p);
		abssub(xt.num.data(), xt.num.size(), q.num.data(), q.num.size());
		while (xt.num.size() > 1 && xt.num.back() == 0) { xt.num.pop_back(); }
		l += p;
		k += p;
		p = xt.num.size() - 1;
		while (!kplan.empty() && k >= kplan.back()) { kplan.pop_back(); }
	}
	xt.addsmall(-1);
	if (floor&&(xt * xt).absbigger(*this, 0)) { xt.addsmall(-1); }
	return xt;
}

integer integer::root(int m,bool floor,int zero)const//zero是浮点补丁
{
	if (!num.back() || m == 1) { return *this; }
	if (m < 1||((m&1)==0&&sign==-1)||(zero&&m==2)) { std::cout << "root"; exit(0); }
	if (m == 2) { return isqrt(floor); }
	int ns = num.size()+zero;
	int f = ns / m + (ns % m != 0);
	const bool useyk = f > 100;
	int k = 1 + (f > 1), need = k + 1, b, b1;
	//view h(*this, m * (f - k));
	//double guess = h.len == 1 ? log(h.ptr[0]) : log(((ll)h.ptr[h.len - 1] << Blen) + h.ptr[h.len - 2]) + (h.len - 2) * log(Base);
	int  hlen=ns-m*(f-k);
	double guess = hlen == 1 ? log(num.back()) : log(((ll)num[num.size()- 1] << Blen) + num[num.size() - 2]) + (hlen - 2) * log(Base);
	integer target = exp(guess / m);
	integer xt = integer(1).shift(k), q = integer(1).shift(k - 1);
	int gap = f > 1 ? Base / 10 : 1;
	ll rate = (ll)1 << 50;
	if (q.absbigger(target, 1)) { target = q + 1; rate >>= 10; }
	else if (target.absbigger(xt, 1)) { target = xt - 1; rate >>= 10; }
	int leftshift;
	//std::cout<<"a\n";
	while (xt.absbigger(q + gap, 0))
	{
		integer mid = (xt + q) / 2;
		if (rate) { mid = (target * (rate - 1) + mid) / rate;}
		b = 0;
		integer tmp = shiftpow(mid, m, need, b);
		leftshift = m * (f - k) + b;
		if (leftshift<zero||(leftshift-zero < ns && !tmp.absbigger(view(*this, leftshift-zero), 0)))
		//由于tmp.size<2*need,而*this在浮点大小root处大小3*l,所以lefshift<zero就一定小
		{
			q = mid;
			if (rate&&q.absbigger(target, 1)) { target = q + 1; rate >>= 10; }
		}
		else
		{
			xt = mid;
			if (rate&&target.absbigger(xt, 1)) { target = xt - 1; rate >>= 10; }
		}
	}
	if (f > 1)
	{
		//std::cout<<"b\n";
		do
		{
			b = 0;
			integer tmp = shiftpow(xt, m - 1, need, b);
			b1 = 0;
			integer tmp1 = shiftmul(xt, tmp, need, b1).shift(b1);
			leftshift = m * (f - k) + b;
			if (ns > leftshift-zero) {
				//std::cout<<f<<" "<<tmp1.num.size()<<"\n";
				view nview(*this, leftshift-zero);
				abssub(tmp1.num.data(), tmp1.num.size(), nview.ptr, nview.len);
				while (tmp1.num.size() > 1 && tmp1.num.back() == 0) { tmp1.num.pop_back(); }
			}
			q = tmp1 / (tmp * m);
		} while (q.sign > 0 && q.num.back() && (xt = xt - q).num.back());
		int p = k - 1;
		integer yt; int l, yk = 0;
		std::vector<int>kplan;
		int k1 = f; k1 += ((k1 & 1) == 1);
		while (k1 > 5) { k1 = k1 / 2 + 1; k1 += ((k1 & 1) == 1); kplan.push_back(k1); }
		bool again=f<10;
		//std::cout<<"c\n";
		while (k < f||again)
		{
			again&=k<f;
			if (p + k > f) { p = f - k;}
			else if (!kplan.empty() && kplan.back() > k && kplan.back() - k < p) { p = kplan.back() - k; kplan.pop_back(); }
			need = xt.num.size() + p + 2, b = 0, b1 = 0;
			integer tmp = shiftpow(xt, m - 1, need, b);
			integer tmp1 = shiftmul(xt, tmp, need, b1).shift(b1 + p);//tmp1最大,决定*this折叠0边界
			leftshift = m * (f - k) + b - p;
			if (ns > leftshift-zero) {
				//std::cout<<f<<" "<<tmp1.num.size()<<"\n";
				view nview(*this, leftshift-zero);
				int check = abssub(tmp1.num.data(), tmp1.num.size(), nview.ptr, nview.len);
				if (check < 0) { std::cout << "root"; exit(0); }
				while (tmp1.num.size() > 1 && tmp1.num.back() == 0) { tmp1.num.pop_back(); }
			}
			bool flag = useyk;
			if (useyk)
			{
				int diff = tmp.num.size() - xt.num.size();
				view tmpview(tmp, diff);
				if (yt.num.empty())
				{
					l = 2 * tmpview.len;
					integer r;
					yt = divide(integer(1).shift(l), tmpview, r);
				}
				else
				{
					yk = tmpview.len - (int)yt.num.size();
					need = yt.num.size() + yk + 2, b = 0;
					integer ytmp = shiftmul(yt, yt, need, b);
					ytmp = shiftmul(tmpview, ytmp, need, b);
					view ytmpview(ytmp, -(yk - l + b));
					yt = (yt + yt).shift(yk);
					abssub(yt.num.data(), yt.num.size(), ytmpview.ptr, ytmpview.len);
					while (yt.num.back() == 0 && yt.num.size() > 1) { yt.num.pop_back(); }
					l += yk;
				}
				need = tmp1.num.size() - tmp.num.size() + 1, b = 0;
				if (need > 7)
				{
					q = shiftmul(tmp1, yt, need, b);
					view qview(q, l + diff - b);
					integer prod = karamul(qview, tmp);
					tmp1.sign = abssub(tmp1.num.data(), tmp1.num.size(), prod.num.data(), prod.num.size());
					while (tmp1.num.back() == 0 && tmp1.num.size() > 1) { tmp1.num.pop_back(); }
					integer q1= divide(tmp1, tmp, tmp1);
					q = addorsub(qview.ptr, qview.len, qview.sign, q1.num.data(), q1.num.size(), q1.sign, 1);
					if (tmp1.sign == -1 && tmp1.num.back()) { q.addsmall(-1); }
					q = q / m;
				}
				else { flag = 0; }
			}
			if (!flag)
			{
				q = divide(tmp1,tmp,tmp1) / m;
			}
			xt = xt.shift(p);
			abssub(xt.num.data(), xt.num.size(), q.num.data(), q.num.size());
			while (xt.num.size() > 1 && xt.num.back() == 0) { xt.num.pop_back(); }
			l += p;
			k += p;
			p = xt.num.size() - 2;
			if(p>2&&p<6&&zero){p--;}//浮点版本的等价again机制
			while (!kplan.empty() && k >= kplan.back()) { kplan.pop_back(); }
		}
	}
	xt.addsmall(-1);
	need = num.size() << 1;
	if (floor&&shiftpow(xt, m, need, b).absbigger(*this, 0)) { xt.addsmall(-1); }
	xt.sign=sign;
	return xt;
}

integer integer::shiftpow(const integer& x, int n, int need, int& b)
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
template <bool WantQ>
integer integer::div_newton(const view& a, const view& b, integer& r,int bzero,bool eps)//bzero,eps是浮点补丁
{
	int la = a.len, lb = b.len+bzero;
	int n =std::min(24.0, la*1.0 / lb+0.2);
	int l = (la - lb) / n + 3,l0=(la - lb) / (n+1) + 3;
	while (l >= 512) { l >>= 1;l0>>=1;}
	if(l>255&&l0<256){n++;}
	l = (la - lb) / n + lb + 2;
	integer xt=reciprocal(b,l-bzero);xt.addsmall(-xt.sign);
	if(r.num.data()!=a.ptr||r.num.size()!=a.len){r=a;}
	integer q; 
	if constexpr (WantQ)
	{
		q.num.assign(la - lb + 1, 0);
	}
    view bview;double*db;
    if(n>1){trycatch();bview.sign=0;}
	for (int i = 0; i < n && r.num.size() > lb + 4; i++) {
		int lr = r.num.size(), start = lr - l;
		if (start < 0) { start = 0; }
		view cutview(r.num.data(), start, lr - 1, r.sign);
		int need = cutview.len - lb + 2, shift = 0;
        if(lastview.sign>0&&need<lastview.len&&need*1.3>lastview.len){need=lastview.len;}
		integer p = shiftmul(cutview,xt, need, shift);
		view pview(p, l - shift);
		if(start==0&&eps)
		{
			q.shiftadd(pview,0);
			if(n>1){endcatch();if(bview.sign>0){delete[]db;}}
			while (q.num.size() > 1 && q.num.back() == 0) { q.num.pop_back(); }
			q.sign = a.sign * b.sign;
			return q;
		}
        view lastviewtmp;double*lastytmp;
        lastviewtmp=lastview,lastytmp=lasty;
        lastview=bview,lasty=db;
		integer now = karamul(pview,b);
        bview=lastview,db=lasty;
        lastview=lastviewtmp,lasty=lastytmp;
		if (abssub(r.num.data() + start+bzero, std::min(l, lr)-bzero, now.num.data(), now.num.size()) < 0) { std::cout << "div_newton1"; exit(0); }
		if constexpr (WantQ)
		{
			q.shiftadd(pview, start);
		}
		while (!r.num.back() && r.num.size() > 1) { r.num.pop_back(); }
	}
    if(n>1){endcatch();if(bview.sign>0){delete[]db;}}
	if (r.num.size() > lb + 4) { std::cout << "div_newton2"; exit(0); }
	if constexpr (WantQ)
	{
		if(bzero){r.div2pow(bzero*Blen);}
		q.shiftadd(div_native<1>(r, b, r), 0);
		while (q.num.size() > 1 && q.num.back() == 0) { q.num.pop_back(); }
		q.sign = a.sign * b.sign;
	}
	else{div_native<0>(r, b, r);}
	r.sign = a.sign;
	return q;
}

integer integer::div_bz(const view& a, const view& b, integer& r)
{
	int m=1,back=b.ptr[b.len-1];
	while(Base>(back<<m)){m++;}m--;
	bool needcopy=r.num.data()!=a.ptr||r.num.size()!=a.len;
	r.num.reserve(a.len+1);
	if(needcopy){r=a;}
	r.mul2pow(m);
	if(r.num.size()==a.len){r.num.push_back(0);}
	integer bb=b;bb.mul2pow(m);
	int bblen=bb.num.size(),bblen2=bblen<<1;
	integer q;q.num.resize(r.num.size()-bblen+1);
	int i=r.num.size()-bblen2;
	while(i>-1)
	{
		div_4n_2n(r.num.data()+i,bblen2,bb.num.data(),bblen,q,i);
		i-=bblen;
		int j=i;
		while(j+bblen2>0&&r.num[j+bblen2-1]==0){j--;}
		if(i!=j){i=j+(r.num[j+bblen2-1]>=bb.num.back());}
	}
	i+=bblen2;
	if(i){div_4n_2n(r.num.data(),i,bb.num.data(),bblen,q,0);}
	r.num.resize(bblen);
	while(r.num.size()>1&&r.num.back()==0){r.num.pop_back();}
	while(q.num.size()>1&&q.num.back()==0){q.num.pop_back();}
	r.div2pow(m);
	r.sign=a.sign;
	q.sign=a.sign*b.sign;
	return q;
}

integer integer::divide(const view& a, const view& b, integer& r, bool wantq)
{
    int la = a.len, lb = b.len;
    if (la < lb) { r = a; return 0; }
    float f = (float)la / lb;
    if (la > 200 && lb > 50 && (la - lb > 100 || f > 3))
    {
        if (lb < 400 - f * 30 && f < 2.6) { return div_bz(a, b, r); }
        if (wantq) { return div_newton<true>(a, b, r); }
        else       { return div_newton<false>(a, b, r); }
    }
    if (wantq) { return div_native<true>(a, b, r); }
    else       { return div_native<false>(a, b, r); }
}

integer integer::mod_positive(const integer& b)const
{
	integer r;
	divide(*this, b, r,0);
	if (r.sign < 0 && r.num.back())
	{
		r.num.resize(b.num.size());
		abssub(r.num.data(), r.num.size(), b.num.data(), b.num.size());
		while (r.num.back() == 0 && r.num.size() > 1) { r.num.pop_back(); }
		r.sign = 1;
	}
	return r;
}

int integer::ctz()const
{
	if(!num.back()){std::cout<<"ctz";exit(0);}
	int i = 0;
	for (; i < num.size() - 1 && num[i] == 0; i++);
	int h = num[i]; i *= Blen;
	while ((h & 1) == 0) { h >>= 1; i++; }
	return i;
}

void integer::mod2pow(int n)
{
	int end = n / Blen; n %= Blen;
	if (end >= num.size()) { return; }
	num[end] &= (1 << n) - 1;
	num.resize(end + 1);
	while(num.size()>1&&!num.back()){num.pop_back();}
}

void integer::div2pow(int n)
{
	int start = n / Blen, i = 0; n %= Blen;
	int m = Blen - n, mask = (1 << n) - 1;
	for (; start + 1 < num.size(); i++)
	{
		int low = num[start] >> n, high = num[++start] & mask;
		num[i] = low | (high << m);
	}
	if (start >= num.size()) { num[i] = 0; }
	else
	{
		int h = num[start] >> n;
		if (h) { num[i] = h; }
		else { num[i] = 0; i -= i > 0; }
	}
	num.resize(i + 1);
}

void integer::mul2pow(int n)
{
	if (!num.back()) { return; }
	int start = n / Blen; n %= Blen;
	int m = Blen - n, mask = (1 << m) - 1, ns = num.size();
	num.resize(start + ns + 1);
	for (int i = ns; i > 0; i--)
	{
		int high = num[i] & mask, low = num[i - 1] >> m;
		num[start + i] = low | (high << n);
	}
	num[start] = (num[0] & mask) << n;
	for (int i = 0; i < start; i++) { num[i] = 0; }
	if (num.back() == 0) { num.pop_back(); }
}

integer integer::shift(int n)const
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

void integer::addsmall(int n)
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

void integer::mul_core(std::vector<integer>& level, int& zeros)
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
#if  shulun>0
// ============ 类外自由函数（前半） ============

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

integer C_core(int e, int s, int  e_s, int& zeros)
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
	for (int i = 0, j = 0; i < u.size() - 1; j = u[i], i += 2)
	{
		l.push_back(power(mul_lambda(j, u[i]), u[i + 1]));
	}
	integer::mul_core(l, s);
	l[0].mul2pow(b%Blen);
	zeros += b / Blen;
	return l[0];
}

integer mul(std::vector<integer>& level)
{
	if(level.empty()){return 1;}
	int zeros = 0; integer::mul_core(level, zeros); return level.back().shift(zeros);
}

integer A_core(int e, int s, int& zeros)
{
	if (s > 10000 && (e >> 2) < s) { return C_core(e, e - s, 0, zeros); }
	std::vector<integer> level; level.reserve(s / 2 + 10);
	for (s = e - s + 1; s < e; s += 2)
	{
		level.push_back((ll)s * (s + 1));
	}
	if (s == e) { level.push_back(e); }
	integer::mul_core(level, zeros);
	return level.back();
}

integer A(int e, int s)
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

integer C(int e, int s)
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

void gcdshift1_inplace(integer&x,integer&y,ll A,ll B,ll C,ll D)
{
	int i = 0; ll kx = 0, ky = 0;
	for (; i < y.num.size(); i++) {
		int xi = x.num[i], yi = y.num[i];
		kx += A * xi + C * yi;
		ky += B * xi + D * yi;
		x.num[i] = kx & Bmask; kx >>= Blen;
		y.num[i] = ky & Bmask; ky >>= Blen;
	}
	for (; i < x.num.size(); i++) {
		kx += A * x.num[i];
		x.num[i] = kx & Bmask; kx >>= Blen;
	}
	while (x.num.size() > 1 && x.num.back() == 0) x.num.pop_back();
	while (y.num.size() > 1 && y.num.back() == 0) y.num.pop_back();
}

void gethigh(integer&x,integer&y,ll&xh,ll&yh)
{
	xh = x.num.back(), yh = y.num.back();
	if (x.num.size() > 1)
	{
		xh = (xh << Blen) + x.num[x.num.size() - 2];
		if (x.num.size() == y.num.size())
		{
			yh = (yh << Blen) + y.num[y.num.size() - 2];
		}
	}
}

bool lehmer(integer& x, integer& y, integer* a, integer* b, integer* c, integer* d)
{
	if (a) { a->sign *= x.sign; c->sign *= x.sign; b->sign *= y.sign; d->sign *= y.sign;}
	x.sign = y.sign = 1;
	ll xh,yh;gethigh(x,y,xh,yh);
	ll A = 1, B = 0, C = 0, D = 1, q, u, v, check;
	while ((u = yh + B) && (v = yh + D) && (check = xh + C - (q = (xh + A) / u) * v) > -1 && check < v)
	{
		u = A, A = B, B = u - B * q;
		u = C, C = D, D = u - C * q;
		u = xh, xh = yh, yh = u - yh * q;
	}
	if (B)
	{
		gcdshift1_inplace(x,y,A,B,C,D);
		auto shift2 = [](integer& a, integer& b, ll A, ll B, ll C, ll D) {
			int i = 0; ll ka = 0, kb = 0;
			if (!b.num.back()) b.sign = -a.sign;
			for (; i < a.num.size(); i++) {
				int ai = a.num[i], bi = b.num[i];
				ka += llabs(A * ai - C * bi);
				kb += llabs(B * ai - D * bi);
				a.num[i] = ka & Bmask; ka >>= Blen;
				b.num[i] = kb & Bmask; kb >>= Blen;
			}
			for (; i < b.num.size(); i++) {
				ka += llabs(C * b.num[i]);
				kb += llabs(D * b.num[i]);
				a.num.push_back(ka & Bmask); ka >>= Blen;
				b.num[i] = kb & Bmask; kb >>= Blen;
			}
			while (ka) { a.num.push_back(ka & Bmask); ka >>= Blen; }
			while (kb) { b.num.push_back(kb & Bmask); kb >>= Blen; }
			b.sign *= (D > 0 ? 1 : -1);
			a.sign = -b.sign;
		};
		if (a) shift2(*a, *b, A, B, C, D), shift2(*c, *d, A, B, C, D);
	}
	return B;
}

int jacobi(const integer& a, const integer& b)
{
	if (b.num[0] % 2 == 0 || b.sign == -1) { std::cout << "jacobi"; exit(0); }
	integer y(a.mod_positive(b)), x(b);
	int m = 1;
	while (y.num.back())
	{
		if (y.num[0] % 2 == 0)
		{
			int k = y.ctz();
			y.div2pow(k);
			int j = x.num[0] % 8;
			int s = (j == 1 || j == 7) ? 1 : -1;
			if (s == -1 && k % 2) { m = -m; }
		}
		if (y.num.size() == 1 && y.num[0] == 1) { return m; }
		int xl=x.num[0]&3,yl=y.num[0]&3;
		if(x.num.size()<y.num.size()+2)
		{
			ll xh,yh;
			gethigh(x,y,xh,yh);
			ll A = 1, B = 0, C = 0, D = 1, q, u, v, check;
			while ((u = yh + B) && (v = yh + D) && (check = xh + C - (q = (xh + A) / u) * v) > -1 &&check < v)
			{
				bool odd=q&1;
				if(odd){if(check<<6<yh){break;}q++;}
				if(xl==3&&yl==3){m=-m;}
				u = A, A = B, B = u - B * q;
				u = C, C = D, D = u - C * q;
				u = xh, xh = yh, yh = u - yh * q;
				u = xl, xl = yl, yl = (u - yl * (q&3))&3;
				if(odd)
				{
					B=-B,D=-D,yh=-yh,yl=4-yl;
					if(xl==3){m=-m;}
				}
			}
			if(B){gcdshift1_inplace(x,y,A,B,C,D);continue;}
		}
		if(xl==3&&yl==3){m=-m;}
		x %=y;
		std::swap(x, y);
	}
	return b.num.size() == 1 && b.num[0] == 1;
}

void hgcd(const integer::view& x, const integer::view& y, integer& a, integer& b, integer& c, integer& d)
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
				integer q = integer::divide(xn, yn, xn);
				gcdshift(a, b, q);
				gcdshift(c, d, q);
				std::swap(xn, yn);
			}
		}return;
	}
	else if (y.len - 1 < m) { return; }
	integer A(1), B(0), C(0), D(1), xn(x), yn(y);
	hgcd(integer::view(xn.num.data(), m, xn.num.size() - 1, xn.sign), integer::view(yn.num.data(), m, yn.num.size() - 1, yn.sign), A, B, C, D);
	gcdshift1(xn, yn, A, B, C, D);
	gcdshift2(a, b, c, d, A, B, C, D);
	if (!yn.num.back() || yn.num.size() - 1 < m) { return; }
	A = 0, B = C = 1, D = integer::divide(xn, yn, xn); D.sign = -D.sign;
	if (yn.num.size() >= 2 * m) { std::cout << "hgcd"; exit(0); }
	if (xn.num.back() && xn.num.size() > m)
	{
		int k = 2 * m - yn.num.size();
		hgcd(integer::view(yn.num.data(), k, yn.num.size() - 1, yn.sign), integer::view(xn.num.data(), k, xn.num.size() - 1, xn.sign), A, B, C, D);
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
			if (!(x.num.size() < y.num.size() + 2 && lehmer(x, y, 0, 0, 0, 0)))
			{
				integer::divide(x, y, x);
				std::swap(x, y);
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
	const bool big = y.num.size() > 900;
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
			integer q=integer::divide(x,y,x);
			gcdshift(a, b, q);
			gcdshift(c, d, q);
			std::swap(x, y);
		}
	}
	if (y.num.back())
	{
		integer A(1), B(0), C(0), D(1);
		while (y.num.back())
		{
			if (!(x.num.size() < y.num.size() + 2 && lehmer(x, y, &A, &B, &C, &D)))
			{
				integer q=integer::divide(x,y,x);
				gcdshift(A, B, q);
				gcdshift(C, D, q);
				std::swap(x, y);
			}
		}
		if (big) { a = a * A + b * C, c = c * A + d * C; }
		else { std::swap(a, A), std::swap(c, C); }
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

integer inv(const integer& a, const integer& p)
{
	integer s, t, g = euclid(p, a, s, t);
	if (!(g.num.size() == 1 && g.num[0] == 1) || !t.num.back()) { return 0; }
	return t.mod_positive(p);
}

integer invbasepower(const integer&a,int m)//a>0的奇数
{
	std::vector<int>mplan;
	mplan.push_back(m);
	while(m>1)
	{
		m=m/2+(m&1);
		mplan.push_back(m);
	}mplan.pop_back();
	unsigned int x1=1,a1=a.num[0];
	for(int i=1;i<32;i<<=1)
	{
		x1=x1*(2-a1*x1);
	}
	integer  x=x1;x.num.resize(1);
	while(!mplan.empty())
	{
		m=mplan.back();mplan.pop_back();
		bool odd=m&1;m+=odd;
		integer tmp=integer::karamul(integer::view(a.num.data(), 0,std::min(m,(int)a.num.size())-1,1),x,0,m);
		while(tmp.num.size()>1&&tmp.num.back()==0){tmp.num.pop_back();}
		tmp.addsmall(-2);
		x=integer::karamul(x,integer::view(tmp.num.data(), 0,std::min(m,(int)tmp.num.size())-1, tmp.sign),0,m-odd);
		x.num.resize(std::min(m-odd,(int)x.num.size()));
		x.sign*=-1;
	}
	return x;
}
// ============ mont / mod 类成员函数 ============

integer::view mont::fastmod(const integer& x)
{
	return  integer::view(x.num.data(), 0, std::min((int)x.num.size() - 1, m - 1), x.sign);
}
void mont::init()
{
	m = p.num.size();
	integer n = integer(1).shift(m);
	a=invbasepower(p,m);
	a.sign*=-1;
	if(a.sign<0){a=a+n;}
	f = n % p;
	b = f * f % p;
}

bool mont::karamulchoosemultiply(integer::view a,integer::view b)
{
	if (a.len < b.len) { std::swap(a, b); }
	bool canfft=a.len + b.len <= halflen;
	if (b.len < 80 || (b.len < 160 &&((!canfft)||a.len * b.len < 30000))) { return 1; }
	return 0;
}

integer mont::out(const integer& x)
{
	if(!x.num.back()){return 0;}
	integer tmp = integer::karamul(fastmod(x), a,0,m);
	integer::view tmpview=fastmod(tmp);
	integer u=integer::karamul(p, tmpview,std::max(0,m-5));
	int carry=0;
	if(m>5&&karamulchoosemultiply(p,tmpview))
	{
		bool add1=0;
		for(int i=0,end=std::min((int)x.num.size(),m-1);i<end;i++){add1|=x.num[i]!=0;}
		int shouldbe=(Base-add1-(x.num.size()<m?0:x.num[m-1]))&Bmask;
		carry=u.num[m-1]>shouldbe;
		//u(m-1)处来自低位量级eps=m*base^2/base^4=0,但是考虑全Bmask影响无法忽略eps=1
		//u(m-1)+carry=shouldbe mod base,由eps=1得carry=0或者1 
		//仅需修正u(m-1)=bmask,shouldbe=0情形,此时向上产生一个进位,判据u(m-1)>shouldbe
		//其余情况不会向上进位
		if(x.num.size()>=m){add1|=x.num[m-1]!=0;}
		carry+=add1;
		int end=m;
		if(u.num.size()<=m){u.num.resize(1);end=1;}
		for(int i=0;i<end;i++){u.num[i]=0;}
	}
	integer t=x+u;
	t.div2pow(m*Blen);
	t.addsmall(carry);
	if (!p.absbigger(t, 0))
	{
		integer::abssub(t.num.data(), t.num.size(), p.num.data(), p.num.size());
		while (t.num.size() > 1 && t.num.back()== 0) { t.num.pop_back(); }
	}
	return t;
}

integer mont::in(const integer& x)
{
	return out((x * b).mod_positive(p));
}

integer mont::add(const integer& u, const integer& v)
{
	integer ans(u + v);
	if (ans.absbigger(p, 1))
	{
		integer::abssub(ans.num.data(), ans.num.size(), p.num.data(), p.num.size());
		while (ans.num.size() > 1 && ans.num.back() == 0) { ans.num.pop_back(); }
	}
	return ans;
}

integer mont::sub(const integer& u, const integer& v)
{
	integer ans;
	if (u.absbigger(v, 1)) { ans = u; }
	else { ans = u + p; }
	integer::abssub(ans.num.data(), ans.num.size(), v.num.data(), v.num.size());
	while (ans.num.size() > 1 && ans.num.back() == 0) { ans.num.pop_back(); }
	return ans;
}

integer mont::pow_binary(const integer& x, const integer& y)
{
	const std::vector<int>& l = y.num;
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

integer mont::pow_exponent(const integer& u, const std::vector<int>& l, int w)
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

// ============ mod 类成员函数 ============

integer mod::pow_exponent_m(const integer& u, const std::vector<int>& l, int w)
{
	if (l.empty()) { return 1; }
	std::vector<integer>b; b.push_back(u);
	integer a = integer::karamul(u,u,0,m); a.mod2pow(m);
	w = (1 << (w - 1)) - 1;
	for (int i = 0; i < w; i++)
	{
		integer tmp = integer::karamul(a,b.back(),0,m);
		tmp.mod2pow(m);
		b.push_back(tmp);
	}
	const int ls = l.size();
	a = b[l[ls - 2]];
	for (int j = 0; j < l[ls - 1]; j++) { a = a * a; a.mod2pow(m); }
	for (int i = ls - 4; i > -1; i -= 2)
	{
		a = integer::karamul(a,b[l[i]],0,m);
		a.mod2pow(m);
		for (int j = 0; j < l[i + 1]; j++)
		{
			a = integer::karamul(a,a,0,m);
			a.mod2pow(m);
		}
	}
	return a;
}

mod::mod(const integer& x)
{
	if (x.num[0] & 1) { a1 = 1, m = 0; a.p = x; a.init(); return; }
	m = x.ctz();
	integer y = x; y.div2pow(m);
	if (a1 = (y.num.size() > 1 || y.num.back() > 1))
	{
		a.p = std::move(y); a.init();
		key=invbasepower(a.p,m/Blen+(m%Blen!=0));
		key.mod2pow(m);
	}
}

pairs mod::in(const integer& x)
{
	pairs ans;
	if (a1)
	{
		ans.f = a.in(x);
	}
	if (m)
	{
		ans.s = integer::view(x.num.data(), 0, std::min((int)x.num.size() - 1, m / Blen), x.sign);
		ans.s.mod2pow(m);
	}
	return ans;
}

pairs mod::mul(const pairs& u, const pairs& v)
{
	pairs ans;
	if (a1)
	{
		ans.f = a.out(u.f * v.f);
	}
	if (m)
	{
		ans.s =integer::karamul(u.s,v.s,0,m);
		ans.s.mod2pow(m);
	}
	return ans;
}

pairs mod::add(const pairs& u, const pairs& v)
{
	pairs ans;
	if (a1)
	{
		ans.f = a.add(u.f, v.f);
	}
	if (m)
	{
		ans.s = u.s + v.s;
		ans.s.mod2pow(m);
	}
	return ans;
}

pairs mod::sub(const pairs& u, const pairs& v)
{
	pairs ans;
	if (a1)
	{
		ans.f = a.sub(u.f, v.f);
	}
	if (m)
	{
		ans.s = u.s - v.s;
		ans.s.mod2pow(m);
	}
	return ans;
}

integer mod::out(const pairs& u)
{
	if (!a1) { return u.s; }
	integer x = a.out(u.f);
	if (!m) { return x; }
	integer tmp =integer::karamul(u.s - x,key,0,m);
	tmp.mod2pow(m);
	return x + a.p * tmp;
}

pairs mod::one()
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

pairs mod::zero()
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

pairs mod::pow_exponent(const pairs& x, const std::vector<int>& l, int w)
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

// ============ exponent / getw / power(x,n,p) ============

std::vector<int> exponent(const integer& a, int w)
{
	std::vector<int>l;
	int step = 0, now = 0, cnt = 0;
	for (int i = 0; i < a.num.size(); i++)
	{
		int m = a.num[i]; bool last = (i + 1 == a.num.size());
		for (int j = 0; j < Blen && (!last || m); j++, m >>= 1)
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

int getw(const  integer& n)
{
	return std::min((int)n.num.size(), 4) + 1;
}

integer power(const integer& x, const integer& n, const integer& p)
{
	if (n.num.size() == 1)
	{
		int m = n.num.back();
		integer a(x.mod_positive(p)), ans = m & 1 ? a : 1;
		while (m >>= 1)
		{
			a = a * a;
			a %=p;
			if (m & 1) { ans = ans * a ;ans%=p; }
		}
		return ans;
	}
	mod q(p); int w = getw(n);
	integer ans = q.out(q.pow_exponent(q.in(x), exponent(n, w), w));
	return ans.sign < 0 && ans.num.back() ? ans + p : ans;
}
// ============ shulun 段：素性 / QS ============
#endif
#if shulun>1

bool inside(int k, const std::initializer_list<int>& t) {
	for (int v : t) {
		if (v == k) { return 1; }
	}
	return 0;
}

bool miller(const integer& n, const std::initializer_list<int> t)
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

bool lucas(const integer& x)
{
	if (x.num.size() == 1 && x.num[0] == 2) { return 1; }
	integer u;
	if (x.num[0] % 2 == 0 ||
		((inside(x.num[0] % 128, { 1, 9, 17, 25, 33, 41, 49, 57, 65, 73, 81, 89, 97, 105, 113, 121 }) &&
		  ((u = x.isqrt(0)) * u).num == x.num))) {
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

integer shanks_sqrt(const integer& a, const integer& p)
{
	int n = 2;
	while (jacobi(n, p) != -1) { n++; }
	integer t = p - 1;
	int e = t.ctz(); t.div2pow(e);
	mont q; q.p = p; q.init();
	int w = getw(t);
	integer y = q.pow_exponent(q.in(n), exponent(t, w), w), a0 = q.in(a),
	x = q.pow_exponent(a0, exponent((t - 1) / 2, w), w), b = q.out(a0 * q.out(x * x));
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

int inv_fsp(int q, int p)
{
	int h=p-2,a = h & 1 ? q : 1;
	while (h >>= 1)
	{
		q = q * q % p;
		if (h & 1) { a = q * a % p; }
	}
	return a;
}

/* ========================= qsinit ========================= */

void qsinit(integer&d,const integer&n,std::vector<int>&prime,std::vector<int>&root,std::vector<uc>&logp,uc&bound,mont&q,int&m)
{
	d = n.isqrt(0);
	int B = pow(d.getlog() - 2, 1.67 + 0.106 * n.num.size());
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
	m = 1 << 15;
	for (int i = 1; i < prime.size(); i++)
	{
		int p = prime[i];
		logp.push_back(log(p)*3);
		root.push_back(shanks_sqrt(n, p).num[0]);
	}
	d = d * 14142 / (integer(m) * 10000);
	bound = (d.getlog() + 2 * log(m) - 13)*3;
	q.p = n; q.init();
}

/* ========================= qsgeta ========================= */

bool qsgeta(std::vector<int>& index, std::vector<int>& prime,std::vector<int>&inva,
			std::vector<int>&root,int&start,int&gap,std::vector<integer>&crtbase,integer&a,integer&d,integer&b)
{
	if (start > std::min((int)prime.size(), 500)) { start = std::min((int)prime.size() / 4, 30); }
	index.clear(); crtbase.clear();
	int i = start;
	integer tmp;tmp.num.push_back(prime[i]);
	do
	{
		index.push_back(i);
		i += gap + rand() % 8;
		if (i + 1 > prime.size())
		{
			gap++;
			if (gap > prime.size() / 10) { gap = 1; }
			start = std::min((int)prime.size() / 4, 30);
			return 0;
		}
		std::swap(a, tmp);
		tmp = a * prime[i];
	} while (d.absbigger(a * prime.back(), 0));
	int j = (d / a).num[0];
	if (prime.size() > 200 && j < 10) { start++; return 0; }
	for (i = 2; i + 1 < prime.size() && prime[i] < j; i++);
	float low = (float)prime[i - 1] / j, high = (float)prime[i] / j;
	if (abs(1 - low) < abs(1 - high)) { i--; }
	for (j = 0; j < index.size() && index[j] != i; j++);
	if (j != index.size()) { start++; return 0; }
	a = a * prime[i];
	index.push_back(i);
	for (i = 0, b = 0; i < index.size(); i++)
	{
		integer p = prime[index[i]], mi = a / p, invi = inv_fsp((mi % p).num[0], p.num[0]), ti = (invi * root[index[i] - 1]) % p * mi % a;
		b = b + ti;
		crtbase.push_back(ti + ti);
	}b = b % a;
	for (int i = 0; i < root.size(); i++)
	{
		inva[i] = inv_fsp((a % prime[i + 1]).num[0], prime[i + 1]);
	}
	start++;
	return 1;
}

/* ========================= qsgetb ========================= */

bool qsgetb(int&cnt, std::vector<int>& index,std::vector<integer>&crtbase,integer&a,integer&b)
{
	cnt++;
	if (cnt + 1 > (1 << (index.size() - 1))) { return 0; }//只要a-b和b中一个
	int u = 1, pos = 0;
	while ((cnt & u) == 0) { u <<= 1; pos++; }
	u = cnt & (u << 1);
	integer& diff = crtbase[pos];
	if (!u) { diff.sign = -1; }
	else { diff.sign = 1; }
	b = (b + diff).mod_positive(a);
	return 1;
}

/* ========================= qsseive ========================= */

void qsseive(std::vector<int>&root,std::vector<int>&inva,std::vector<int>&prime,std::vector<uc>&logp,std::vector<uc>&pos,std::vector<uc>&neg,integer&b,int m)
{
	for (int i = 0; i < root.size(); i++)
	{
		int invai = inva[i]; if (!invai) { continue; }
		int p = prime[i + 1], h = (b % p).num[0];
		int x0 = (-h - root[i]) * invai % p;
		int x1 = (-h + root[i]) * invai % p;
		if (x0 < 0) { x0 += p; }if (x1 < 0) { x1 += p; }
		const uc lp = logp[i];
		for (int j = x0; j < m; j += p) { pos[j] += lp; }
		for (int j = x1; j < m; j += p) { pos[j] += lp; }
		for (int j = p - x0; j < m; j += p) { neg[j] += lp; }
		for (int j = p - x1; j < m; j += p) { neg[j] += lp; }
	}
}
// ============ factor / primepow / euler / mobius / primeroot / order / fib ============

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
	if (n.getlog()>120) { std::cout << "qs too big"; exit(0); }
	//qs知乎学的,10^47之内保证不超过10s
	//11111111111111111111111111122332231111111111111117791这个需要20s,通过euler的参数l调用,不要直接调用factor
	mont q;
	std::vector<int>prime;
	std::vector<int>root;
	std::vector<uc>logp;
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
	std::vector<uc>pcnt;
	int ta = 0, tb = 0, tc = 0;uc bound;
	auto insert_lambda =
	[&q, &prime, &smooth, &flaw, &fs, &n, &ok, &ta, &tb, &tc, &bound, &index, &pcnt](const integer& x, integer y) -> integer {
		int  ps = prime.size();
		equa e; e.mask.assign(ps / 64 + (bool)(ps % 64 != 0), 0);
		int z = 0;
		integer p=0,res;
		for (; z < ps && (y.num.size() > 1 || y.num[0] > 1); z++)
		{
			p.num[0]=prime[z];
			integer q = integer::div_native<1>(y, p, res);
			int cnt = 0;
			while (!res.num.back())
			{
				cnt++;
				std::swap(y, q);
				q = integer::div_native<1>(y, p, res);
			}
			pcnt[z] = cnt;
		}
		if (y.num.size() > 1 || (y.num[0] > (Base >> 4))) { tb++; return 0; }
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
			if (t) { e.y = e.y * power(prime[i], t); }
		}
		e.y%=n;
		if (y.sign < 0) { e.y = e.y + n; }
		if (y.num[0] > 1)
		{
			equa& s = flaw[j];
			if (s.x.sign < 0) { e.x.sign = y.num[0]; flaw[j] = e; fs++; return 0; }
			tc++;
			s.x.sign = 1; e.x = q.out(e.x * s.x); s.x.sign = y.num[0];
			e.y = (e.y * y.num[0]) * s.y;integer::div_native<0>(e.y,n,e.y);
			for (c = 0; c < cend; c += 64)
			{
				unsigned long long ei = e.mask[c / 64], si = s.mask[c / 64], u = ei & si;
				for (r = 0; r < ps - c && r < 64; r++, u >>= 1)
				{
					if (u & 1) { e.y = e.y * prime[c + r];integer::div_native<0>(e.y,n,e.y);}
				}
				e.mask[c / 64] ^= si;
			}
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
			e.y=e.y*s.y;integer::div_native<0>(e.y,n,e.y);
			for (; c < cend; c += 64)
			{
				unsigned long long ei = e.mask[c / 64], si = s.mask[c / 64], u = ei & si;
				for (r = 0; r < ps - c && r < 64; r++, u >>= 1)
				{
					if (u & 1)
					{
						e.y=e.y*prime[c+r];
						integer::div_native<0>(e.y,n,e.y);
					}
				}
				e.mask[c / 64] ^= si;
			}
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
		if (tb > (ta << (level + 1))) { bound ++; }if ((ta << level) > tb) { bound --; }//倾向于错误多些,快
		return 0;
	};
	int m;integer d;
	qsinit(d, n, prime, root, logp, bound, q,m);
	pcnt.resize(prime.size());
	smooth.assign(prime.size(), equa(n.num.size(), prime.size()));
	flaw.assign(prime.size() << 2, equa(n.num.size(), prime.size()));
	ok = 40;
	integer a, b; int start = 30, gap = 1;
	std::vector<integer>crtbase;
	std::vector<int>inva; inva.resize(root.size());
	int cnt = 0;
	std::vector<uc>pos(m,0), neg(m,0);
	while (ok > 0)
	{
		if (!qsgetb(cnt,index,crtbase,a,b)) { cnt = 0; while (!qsgeta(index,prime,inva,root,start,gap,crtbase,a,d,b));}
		integer c = (b * b - n) / a;
		qsseive(root, inva, prime, logp, pos, neg, b, m);
		for (int i = 0; i < m; i++)
		{
			if (pos[i] > bound)
			{
				g = insert_lambda(b + a * i, (a * ((ll)i * i) + b * (i << 1) + c));
				if (check_lambda()) { return g; }
			}pos[i] = 0;
			if (neg[i] > bound)
			{
				g = insert_lambda(b - a * i, (a * ((ll)i * i) - b * (i << 1) + c));
				if (check_lambda()) { return g; }
			}neg[i] = 0;
		}
	}
	std::cout << "factor"; exit(0);//可能是输入质数幂(qs非法输入)，或者极小概率geta产生重复的a导致，或者ok正常耗尽
}

integer primepow(const integer& x)
{
	integer p;
	if (x.num[0] % 2 == 0) { p = 2; }
	else
	{
		int w = getw(x);
		std::vector<int>b = exponent(x, w);
		ll a = 2;
		mont q; q.p = x; q.init();
		do
		{
			p = gcd(q.out(q.pow_exponent(q.in(a), b, w)) - a, x);
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

integer euler(const integer& x, std::vector<integer>* l, bool p)
{
	if (x.num.size() == 1 && x.num[0] < 2) { return 1; }
	integer g = primepow(x);
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

integer euler(const integer&n,std::vector<integer>&l,std::vector<integer>&l0)
{
	integer phi=euler(n,&l);
	for(int i=0;i<l.size();i++)
	{
		integer g=gcd(phi,l[i]);
		if(g.num.size()>1||g.num[0]>1)
		{
			l0.push_back(l[i]);
		}
	}
	for(int i=0;i<l.size();i++)
	{
		std::vector<integer>tmp;
		euler(l[i]-1,&tmp);
		for(auto&e:tmp)
		{
			int j=0;
			for(;j<l0.size()&&l0[j].num!=e.num;j++);
			if(j==l0.size()){l0.push_back(e);}
		}
	}
	return phi;
}

int mobius(const integer& x, bool p)
{
	if (x.num.size() == 1 && x.num[0] == 1) { return 1; }
	integer g = primepow(x);
	if (g.num.back()) { return x.num == g.num ? -1 : 0; }
	integer u = factor(x, p), v = x / u;
	g = gcd(u, v);
	if (g.num.size() > 1 || g.num[0] > 1) { return 0; }
	int a = mobius(u, p);
	return a ? a * mobius(v, p) : 0;
}

int primeroot(const integer& p,const std::vector<integer>*k0)
{
	if (p.num.size() == 1 && p.num[0] < 6) { return p.num[0] % 2 + 1; }
	std::vector<integer>k;
	integer x= p - 1;
	if(k0){k=*k0;}
	else{euler(x, &k);}
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
	l.push_back(k[0]);
	for (int i = 1; i < k.size(); i++)
	{
		l.push_back(k[i] - k[i - 1]);
	}
	mont q; q.p = p; q.init();
	for (int a = 2;; a++)
	{
		integer t = q.in(a), y = q.pow_binary(t, l[0]);
		int s = 0;
		while (y.num != q.f.num)
		{
			if (++s == l.size()) { return a; }
			y = q.out(y * q.pow_binary(t, l[s]));
		}
	}
}

integer order(const integer&a,const integer&n,integer phi,const std::vector<integer>&l0)
{
	integer g=gcd(a,n);
	if(g.num.size()>1||g.num[0]>1){return 0;}
	for(auto&e:l0)
	{
		integer r,next=integer::div_native<1>(phi,e,r);
		do
		{
			integer check=power(a,next,n);
			if(check.num.size()>1||check.num[0]>1){break;}
			std::swap(phi,next);
			next=integer::divide(phi,e,r);
		} while (r.num.back()==0);
	}
	return phi;
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
// ============ indexcalculus / dlsolver + 收尾 ============

/* ========================= indexcalculus ========================= */

bool indexcalculus::smooth(integer tmp,std::vector<int>&right)
{
	right.clear();
	int j = 0, ps = prime.size();
	integer p=0,r;
	for (; j < ps && (tmp.num.size() > 1 || tmp.num[0] > 1); j++)
	{
		p.num[0]=prime[j];
		integer q = integer::div_native<1>(tmp, p, r);
		int cnt = 0;
		while (r.num.back() == 0)
		{
			cnt++;
			std::swap(tmp, q);
			q = integer::div_native<1>(tmp, p, r);
		}
		if (cnt) { right.push_back(j), right.push_back(cnt); }
	}
	return j < ps;
}

void indexcalculus::prepare()
{
	double lg = q.p.getlog();//lnp<42保证7s内
	if(lg<20||lg>55){std::cout<<"indexcalculus reject";exit(0);}
	if(lg>40){
		std::cout << "\nif mingw64 -O2,ic estimated init time :" << 0.1216* exp(0.2427 *lg) / 1000.0 << " s\n";
	}//   Performance models (-O2,empirical):
	int B=(15-lg*0.2)*pow(1.135,lg);
	bound=lg*4-10;
	flag.assign((B + B % 2) / 2, 1);
	for (int i = 1, a; 2 * (a = 2 * i * (i + 1)) <= B - 1; i++)
	{
		if (flag[i])
		{
			for (int j = a; j < flag.size(); j += 2 * i + 1) { flag[j] = 0; }
		}
	}prime.push_back(2);
	for (int i = 1; i < flag.size(); i++) { if (flag[i]) { prime.push_back(2 * i + 1); } }
	flag.resize(prime.size());
	for (int i = 0; i < prime.size(); i++) { flag[i] = inv_fsp((q.p % prime[i]).num[0], prime[i]); }
	logp.resize(prime.size());
	for (int i = 0; i < prime.size(); i++) { logp[i] = log(prime[i])*4; }
	m = 1 <<(12+q.p.num.size());
}


indexcalculus::indexcalculus(const integer& p0,int r0,const integer&cutoff)
//传入质数p0,原根r,系统计算logr(pi)mod ((p0-1)/cutoff)
//p-1的质因子不能太多,否则使用dlsolver
{
	q.p=p0,pp=(q.p-1)/cutoff;init(r0);
}

integer ax_by_pp(const integer&a,const integer&x,const integer&b,const integer&y,std::vector<ll>&ans,const integer&pp)
{
	int la=a.num.size(),lx=x.num.size(),lb=b.num.size(),ly=y.num.size();
	ans.assign(std::max(la+lx,lb+ly),0);
	for(int i=0;i<la;i++)
	{
		for(int j=0;j<lx;j++)
		{
			ans[i+j]+=(ll)a.num[i]*x.num[j];
		}
	}
	integer c;c.sign=a.sign*x.sign;
	bool absadd=c.sign==b.sign*y.sign;
	if(absadd)
	{
		for(int i=0;i<lb;i++)
		{
			for(int j=0;j<ly;j++)
			{
				ans[i+j]+=(ll)b.num[i]*y.num[j];
			}
		}
	}
	else
	{
		for(int i=0;i<lb;i++)
		{
			for(int j=0;j<ly;j++)
			{
				ans[i+j]-=(ll)b.num[i]*y.num[j];
			}
		}
	}
	c.num.reserve(ans.size()+absadd);
	ll carry=0;
	for(int i=0;i<ans.size();i++)
	{
		carry+=ans[i];
		c.num.push_back(carry&Bmask);
		carry>>=Blen;
	}
	while(carry>0)
	{
		c.num.push_back(carry&Bmask);
		carry>>=Blen;
	}
	if(carry<0){c.num.back()-=Base;}
	while(c.num.size()>1&&c.num.back()==0){c.num.pop_back();}
	if(c.num.back()<0)
	{
		int k=1;
		for(int i=0;i<c.num.size();i++)
		{
			k+=Bmask-c.num[i];
			c.num[i]=k&Bmask;
			k>>=Blen;
		}
		while(c.num.size()>1&&c.num.back()==0){c.num.pop_back();}
		c.sign*=-1;
	}
	integer::div_native<0>(c,pp,c);
	return  c;
}

void indexcalculus::init(int r0)
{
	prepare();
	int ps = prime.size();
	struct equa
	{
	public:
		std::vector<integer>h;
		void clear(){for (auto& e : h) { e.num.resize(1),e.num[0]=0;}}
		equa(int ps, const integer& p){h.assign(ps + 1, p);clear();}
	};
	std::vector<equa>matrix;
	matrix.assign(ps, equa(ps, pp));
	std::vector<bool>need(ps, 1); int index = 0;
	std::vector<bool>needb(ps, 1); int indexb = 0, indexbb = 1, indexbbb = ps - 1;
	std::vector<int>left,right;
	int maxtry = 24*q.p.num.size()+20,ok=ps+40;
	equa a(ps, pp);std::vector<ll>ans;
	bool needr=1;
	auto insert_lambda = [&a,&needr,&matrix, &ps, &need, &index, &needb, &indexb, &indexbb, &indexbbb, &maxtry, &left,&right,&ok,&ans](const  integer& pp)->void {
		if(--ok<0){std::cout<<"p-1 too many factors,use dlsolver";exit(0);}
		for (int i = 0; i < right.size(); i += 2)
		{
			a.h[right[i]]=right[i + 1];
		}
		for (int i = 0; i < left.size(); i++)
		{
			a.h[left[i]].addsmall(-1);
		}
		a.h[ps] = needr;needr=0;
		int first = 0;
		while (a.h[first].num.back() == 0 && first < ps) { first++; }
		while (index==ps?first <= indexbbb:(first < ps && !need[first]))
		{
			equa& b = matrix[first];
			if (needb[first])
			{
				integer x, y, g = euclid(a.h[first], b.h[first], x, y);
				if (g.num != b.h[first].num)
				{
					b.h[first] = g; g = gcd(pp, g);
					needb[first] = g.num.size() > 1 || g.num[0] > maxtry;
					for (int i = first + 1; i < ps + 1; i++)
					{
						//integer tmp= (a.h[i] * x + b.h[i] * y);
						//integer::div_native<0>(tmp,pp,tmp);
						//b.h[i] =tmp;
						b.h[i]=ax_by_pp(a.h[i],x,b.h[i],y,ans,pp);
					}
				}
				//if(!needb[first]){std::cout<<first<<" good\n";}
			}
			a.h[first].sign*=-1;
			for (int i = first + 1; i < ps + 1; i++)
			{
				//integer tmp=(a.h[i] * b.h[first] - b.h[i] * a.h[first]);
				//integer::div_native<0>(tmp,pp,tmp);
				//a.h[i] = tmp;
				a.h[i]=ax_by_pp(a.h[i],b.h[first],b.h[i],a.h[first],ans,pp);
			}
			a.h[first].sign*=-1;
			first++;
			while (a.h[first].num.back() == 0 && first < ps) { first++; }
		}
		if (first < ps&& index < ps)
		{
			need[first] = 0;
			integer g = gcd(pp, a.h[first]);
			needb[first] = g.num.size() > 1 || g.num[0] > maxtry;
			//if(!needb[first]){std::cout<<first<<" good\n";}
			for (int i = first; i < ps + 1; i++) { matrix[first].h[i] = a.h[i]; }
			for (; index < ps && !need[index]; index++);
		}
		for (; indexb < index && !needb[indexb]; indexb++);
		for (indexbb = indexb + 1; indexbb < index && !needb[indexbb]; indexbb++);
		for (; indexbbb > indexb + 1 && !needb[indexbbb]; indexbbb--);
		a.clear();
	};
	if(smooth(r0,right)){insert_lambda(pp);needr=0;}
	std::vector<uc>pos(m,0); int start = 0;
	while (indexb < ps)
	{
		//std::cout<<indexb<<" "<<" "<<index<<" "<<ps<<"\n";
		left.clear();
		ll h = prime[indexb]; left.push_back(indexb);
		if (index < ps) { h = h * prime[index]; left.push_back(index); }
		else if (indexbb < index) { h = h * prime[indexbb]; left.push_back(indexbb); }
		integer u = needr?h*r0:h; start = (start + 1) % ps;
		for (int j = 0, i = start; j < 3; j++)
		{
			u = u * prime[i];
			left.push_back(i);
			i = (i + 1 + rand() % 8) % ps;
		}
		for (int i = 0; i < ps; i++)
		{
			int pi = prime[i];
			for (int j = pi - ((ll)flag[i] * (u % pi).num[0] % pi); j < m; j += pi) { pos[j] += logp[i]; }
		}
		for (int i = 1; i < m; i++)
		{
			if (pos[i] > bound&&smooth(u + q.p * i,right)) 
			{ 
				insert_lambda(pp);
				if(indexb==ps){break;}
			}
			pos[i] = 0;
		}
	}
	q.init();
	integer cutoff=(q.p-1)/pp;pr=r0;
	integer r=q.pow_binary(q.in(r0),cutoff);
	dlogp.assign(ps, pp);
	for (int i = ps - 1; i > -1; i--)
	{
		integer tmp = matrix[i].h[ps];
		for (int j = ps - 1; j > i; j--)
		{
			tmp = tmp - dlogp[j] * matrix[i].h[j];
		}
		integer end = q.pow_binary(q.in(prime[i]),cutoff),s,t,g =euclid(pp, matrix[i].h[i],s,t), ppp = pp / g;
		tmp = ((tmp / g) *t).mod_positive(ppp);
		s = q.pow_binary(r, tmp), t = q.out(q.pow_binary(r, ppp - tmp) * s);
		int cnt = 0;
		while (cnt<g.num[0]-1&&s.num != end.num)
		{
			s = q.out(s * t); cnt++;
		}
		dlogp[i] = tmp + ppp * cnt;
		//if(q.pow_binary(r,dlogp[i]).num!=end.num){std::cout<<"bad";exit(0);}
	}
}

integer indexcalculus::dlogr(const integer& x)//不要传0
{
	std::vector<int>right;
	integer x0=x.mod_positive(q.p);
	if(smooth(x0,right))
	{
		integer tmp = 0;
		for (int j = 0; j < right.size(); j += 2)
		{
			tmp = tmp + dlogp[right[j]] * right[j + 1];
		}
		return  tmp.mod_positive(pp);
	}
	int cnt=1;
	std::vector<uc>pos(m, 0);
	while(1)
	{
		x0=(x0*pr)%q.p;
		for (int i = 0; i < prime.size(); i++)
		{
			int pi = prime[i];
			for (int j = pi - ((ll)flag[i] * (x0 % pi).num[0] % pi); j < m; j += pi) { pos[j] += logp[i]; }
		}
		for (int i = 1; i < m; i++)
		{
			if (pos[i] > bound && smooth(x0 + q.p * i,right))
			{
				integer tmp = cnt;tmp.sign=-1;
				for (int j = 0; j < right.size(); j += 2)
				{
					tmp = tmp + dlogp[right[j]] * right[j + 1];
				}
				return  tmp.mod_positive(pp);
			}
			pos[i] = 0;
		}cnt++;
	}
}

int dlsolver::bsgs(const integer&a,const integer&b,int m)
{
	if(b.num==ic.q.f.num){return 0;}
	if(a.num==b.num){return 1;}
	integer c=ic.q.pow_binary(a,m),now=c;
	std::vector<int>t(m);
	for(int j=0;j<m;j++)
	{
		t[j]=now.num[0];now=ic.q.out(now*c);
	}now=b;
	for(int i=0;i<m;i++)
	{
		for(int j=0;j<m;j++)
		{
			if(t[j]==now.num[0]&&ic.q.pow_binary(c,j+1).num==now.num){return (j+1)*m-i;}
		}
		now=ic.q.out(now*a);
	}
	std::cout<<"bsgs";exit(0);
}

integer dlsolver::hellman(const integer&a,const integer&b)
{
	std::vector<int>s;
	integer a0=ic.q.pow_binary(ic.q.in(a),ic.pp),b0=ic.q.pow_binary(ic.q.in(b),ic.pp),o=ord(a0,&s);
	integer u=0,v=1;
	for(int i=0;i<s.size();i+=2)
	{   
		int ei=s[i+1];
		integer  pi=s[i],w=power(pi,ei),unit=ic.q.pow_binary(a0,o/pi);
		integer h=0,basep=1;
		int m=sqrt(s[i])+1;
		for(int j=0;j<ei;j++)
		{
			integer f=ic.q.pow_binary(a0,h);
			f=ic.q.pow_binary(f,ic.q.p-2);
			f=ic.q.out(f*b0);
			int ans=bsgs(unit,ic.q.pow_binary(f,o/power(pi,j+1)),m);
			h=basep*ans+h;
			basep=basep*pi;
		}
		u=u+v*(inv(v,w)*(h-u)%w);
		v=v*w;
	}
	return u.mod_positive(o);
}

dlsolver::dlsolver(integer p)
{
	integer small=p.num.size()>2?2147483647:1000000;
	ic.q.p=p;p.addsmall(-1);
	euler(p,&k);
	for(auto&e:k)
	{
		if(small.absbigger(e,0))
		{
			int cnt=0;
			integer r,q=integer::div_native<1>(p,e,r);
			do
			{
				std::swap(p,q);
				q=integer::div_native<1>(p,e,r);cnt++;
			}while(r.num.back()==0);
			l.push_back(e.toll());l.push_back(cnt);
		}
	}
	ic.pp=p;
	if(p.num.size()>1||p.num[0]>1)
	{
		ic.init(primeroot(ic.q.p,&k));
	}
	else{ic.q.init();}
}

integer dlsolver::ord(const integer&a,std::vector<int>*s)
{
	integer a0=s?a:ic.q.in(a),o=s?(ic.q.p-1)/ic.pp:ic.q.p-1;
	for(int i=0;i<l.size();i+=2)
	{
		int cnt=0;
		integer e=l[i],r,next=integer::div_native<1>(o,e,r);
		do
		{
			if(ic.q.pow_binary(a0,next).num!=ic.q.f.num){break;}
			std::swap(o,next);
			next=integer::divide(o,e,r);cnt++;
		} while (r.num.back()==0);
		if(s&&cnt!=l[i+1]){s->push_back(l[i]),s->push_back(l[i+1]-cnt);}
	}
	return o;
}

integer dlsolver::dlogp(const integer&a,const integer&b)
{
	if(ic.q.pow_binary(ic.q.in(b),ord(a)).num!=ic.q.f.num){return -1;}
	if(ic.prime.empty()){return hellman(a,b);}
	integer ans=ic.dlogr(b)*inv(ic.dlogr(a),ic.pp)%ic.pp;
	if(!l.empty())
	{
		integer cutoff=(ic.q.p-1)/ic.pp;
		return ans+(inv(ic.pp,cutoff)*(hellman(a,b)-ans)).mod_positive(cutoff)*ic.pp;
	}
	return ans;
}
#endif
#if fudian
integer floatb::floatbtmp;
int floatb::minvalidlength=500;
integer::view floatb::getview() const
{
    return integer::view(NUM.num.data(), LOW, NUM.num.size() - 1, NUM.sign);
}

integer::view floatb::leftshiftview(int&zero,int leftshift) const
{
	integer::view a=getview();
	if(leftshift>ZERO){leftshift-=ZERO;zero=0;return integer::view(a,leftshift);}
	else{zero=ZERO-leftshift;leftshift=0;return a;}
}

integer::view floatb::clearview(int&zero,int&exp,int l) const//丢弃指数小于l后的视图
{
	integer::view a=getview();
	int aexp=EXP,azero=ZERO;
	if(aexp>=l){zero=azero,exp=aexp;return a;}
	exp=l;
	if(aexp+azero>=l){zero=azero-l+aexp;return a;}
	if(aexp+azero+a.len>l){zero=0;return integer::view(a,l-aexp-azero);}
	zero=-1;return a;
}

floatb floatb::addorsub(const floatb& a, const floatb& b, bool add,int l)//进一步舍弃结果exp<l的部分
{
	if(l<a.EXP){l=a.EXP;}
	if(l<b.EXP){l=b.EXP;}
	int azero,bzero,aexp,bexp;
    integer::view aview = a.clearview(azero,aexp,l), bview = b.clearview(bzero,bexp,l);
	floatb ans;
	if(azero==-1)
	{
		if(bzero==-1){ans.EXP=l,ans.NUM=0;return ans;}
		ans.EXP=bexp;ans.ZERO=bzero;
		ans.NUM=bview;
		ans.NUM.sign*=add?1:-1;
		return ans;
	}
	if(bzero==-1){ans.NUM=aview;ans.EXP=aexp;ans.ZERO=azero;return ans;}
	int change=1;
	if(aexp<bexp){std::swap(aview,bview);std::swap(aexp,bexp);std::swap(azero,bzero);if(!add){change*=-1;}}
	const bool absadd=(aview.sign==bview.sign)==add;
	int gap=aexp-bexp;
	if(gap>=bview.len+bzero){ans.ZERO=azero;ans.EXP=aexp;ans.NUM=aview;ans.NUM.sign*=change;return ans;}
	if(bzero>gap){bzero-=gap;}
	else{gap-=bzero;bzero=0;bview=integer::view(bview,gap);}
	if(azero+aview.len<bzero+bview.len){std::swap(aview,bview);std::swap(azero,bzero);if(!add){change*=-1;}}
	ans.ZERO=std::min(azero,bzero);
	ans.EXP=aexp;
	ans.NUM.num.reserve(azero-ans.ZERO+aview.len+absadd);
	ans.NUM.num.assign(azero-ans.ZERO,0);
	for(int i=0;i<aview.len;i++){ans.NUM.num.push_back(aview.ptr[i]);}
	ans.NUM.sign=aview.sign*change;
    if(absadd)
	{
		ans.NUM.num.push_back(0);
		ans.NUM.shiftadd(bview,bzero-ans.ZERO);
		if(!ans.NUM.num.back()){ans.NUM.num.pop_back();}
	}
	else 
	{
		if(integer::abssub(ans.NUM.num.data(),ans.NUM.num.size(),bview.ptr,bview.len,bzero-ans.ZERO)<0){ans.NUM.sign*=-1;}
		while(ans.NUM.num.size()>1&&ans.NUM.num.back()==0){ans.NUM.num.pop_back();}
	}
	ans.foldzero();
	return ans;
}

void floatb::quotient(const integer::view& a, const integer::view& b,int l,int azero,int  bzero, integer& ans)
{
    int m = 1, t = b.ptr[b.len - 1];
    while ((t << m) < Base) { m++; }
    m -= 1;
    int m_5 = m + 5;
    ull bb = ((ull)b.ptr[b.len - 1] << m_5)
           + (((ull)b.ptr[b.len - 2] << m) >> (Blen - 5));
    floatbtmp.num.reserve(a.len +azero+1);
	floatbtmp.num.assign(azero,0);
	for(int i=0;i<a.len;i++){floatbtmp.num.push_back(a.ptr[i]);}
    floatbtmp.num.push_back(0);
    int ln,w=l-b.len-bzero-1;
	if(w>-1){std::cout<<"quotient request out of range";exit(0);}
    for (;(t =(ln = floatbtmp.num.size())+w) > 0&&ln>2; floatbtmp.num.pop_back())
    {
        int q = (((((ull)floatbtmp.num.back() << Blen) + floatbtmp.num[ln - 2]) << m_5)
                + (((ull)floatbtmp.num[ln - 3] << m) >> (Blen - 5))) / bb;
        if (!q) { continue; }
        int* start = floatbtmp.num.data();
        ll k = 0;
        int end = ln-1, gap = b.len - end,istart=std::max(0,-gap);
        for (int i =istart; i < end; i++)
        {
            k += (ll)start[i] - (ll)b.ptr[i + gap] * q;
            start[i] = k & Bmask;
            k >>= Blen;
        }
        start[end] += k;
        int top = end;
        while (top > -1 && start[top] == 0) { top--; }
        if (top > -1 && start[top] < 0)
        {
            int carry = 0;
            for (int i = istart; i < end; i++)
            {
                carry += start[i] + b.ptr[i + gap];
                start[i] = carry & Bmask;
                carry >>= Blen;
            }
            q--;
        }
        ans.num[t] = q;
    }
	if(t>-1)
	{
		ll u=((ll)floatbtmp.num.back()<< Blen) + floatbtmp.num[ln-2];
		ans.num[t]=u/b.ptr[b.len-1];
		if(t>0)
		{
			int v=u-(ll)b.ptr[b.len-1]*ans.num[t];
			ans.num[0]=v/b.ptr[b.len-1];
		}
	}
}

integer floatb::div_native(const integer::view& a, const integer::view& b,int l,int azero,int  bzero)
{
	integer ans;
	ans.num.resize(l+1);
    if (b.len == 1)
    {
        ull n = b.ptr[0], now = 0;
		int i=a.len-1,gap=azero+l-bzero,end=-gap-1;
        for (int end1=std::max(-1,end); i >end1; i--)
        {
            now = (now << Blen) + a.ptr[i];
            ans.num[i+gap] = now / n;
            now %= n;
        }
		for (;i >end&&now; i--)
        {
            now<<=Blen;
            ans.num[i+gap] = now / n;
            now %= n;
        }
    }
	else{quotient(a, b,l,azero,bzero,ans);}
    ans.sign = a.sign * b.sign;
    while (ans.num.size() > 1 && ans.num.back() == 0) { ans.num.pop_back(); }
	return ans;
}

integer floatb::divide_core(const integer::view& a, const integer::view& b,int l,int azero,int  bzero)
{
	if(a.len+azero!=b.len+bzero||l>b.len+bzero){std::cout<<"divide_core";exit(0);}
	if(l>50&&b.len>50&&(b.len>bzero?(b.len>700-bzero*550.0/b.len):(b.len+bzero>300)))
	{
		int n=azero+l;
		if(n>bzero)
		{
			n-=bzero;
			floatbtmp.num.reserve(n+a.len);
			floatbtmp.num.assign(n,0);
			for(int i=0;i<a.len;i++){floatbtmp.num.push_back(a.ptr[i]);}
			floatbtmp.sign=a.sign;
			return integer::div_newton(floatbtmp,b,floatbtmp,0,1);
		}
		else
		{
			bzero-=n;
			floatbtmp=a;
			return integer::div_newton(floatbtmp,b,floatbtmp,bzero,1);
		}
	}
	else
	{
		return div_native(a,b,l,azero,bzero);
	}
}

floatb floatb::multiply(const floatb& b,int l) const//尝试结果validlehgth=l
{
	integer::view aview = getview(), bview = b.getview();
	int azero=ZERO,bzero=b.ZERO,la = aview.len+azero, lb = bview.len+bzero;
	if(l>la){l=la;}
	if(l>lb){l=lb;}
	if(l<1){l=1;}
	floatb ans;ans.EXP=EXP+b.EXP;
	int gap=la-l;
	if(gap>0){aview=leftshiftview(azero,gap);ans.EXP+=gap;}
	gap=lb-l;
	if(gap>0){bview=b.leftshiftview(bzero,gap);ans.EXP+=gap;}
	ans.EXP+=azero+bzero;
	ans.NUM = integer::shiftmul(aview,bview,l+ 2, ans.EXP);
    ans.setvalidlength(l);
	ans.foldzero();
	return ans;
}

floatb floatb::divide(const floatb& b,int l) const//尝试NUM层面a*base^l/b,尝试结果validlehgth=l
{
	if(!b.NUM.num.back()){std::cout<<"/0";exit(0);}
	integer::view aview = getview(), bview = b.getview();
	int azero=ZERO,bzero=b.ZERO,la = aview.len+azero, lb = bview.len+bzero;
	if(l>la){l=la;}
	if(l>lb){l=lb;}
	if(l<1){l=1;}
	floatb ans;ans.EXP=EXP-b.EXP-l;
	int gap=la-l;
	if(gap>0){aview=leftshiftview(azero,gap);ans.EXP+=gap;}
	gap=lb-l;
	if(gap>0){bview=b.leftshiftview(bzero,gap);ans.EXP-=gap;}
	ans.NUM=divide_core(aview,bview,l,azero,bzero);
	ans.setvalidlength(l);
	ans.foldzero();
	return ans;
}

floatb floatb::calpower(int n,int l)const//指数容易静默溢出,超大指数floatb乘法也是,不检查
{
	floatb c=integer(NUM.num.back()?1:0);
	int l0=validlength();
	if(l>l0){l=l0;}
	if(l<1){l=1;}
	c.setvalidlength(l);
	if(!n||!NUM.num.back()){return c;}
	floatb a;
	if(n<0){a=c.divide(*this,l);n=-n;}
	else{a=*this;a.setvalidlength(l);}
	floatb ans=n&1?a:c;
	while (n >>= 1)
	{
		a = a.multiply(a,l);
		if (n & 1) { ans = ans.multiply(a,l); }
	}
	ans.foldzero();
	return ans;
}

floatb floatb::calroot(int n,int l)const//validlenth()*n要在int一定余量以内,否则溢出,之前太多int了,改不过来了,范围就小点吧
{
	floatb c=integer(NUM.num.back()?1:0);
	int l0=validlength();
	if(l>l0){l=l0;}
	if(l<1){l=1;}
	c.setvalidlength(l);
	if(!n||!NUM.num.back()){return c;}
	floatb a;
	if(n<0){a=c.divide(*this,l);n=-n;}
	else{a=*this;}
	a.setvalidlength(l);
	if(n==1){return a;}
	ll t=l*n,gap=(t-l)*Blen,h=(n-(gap-(ll)a.EXP*Blen)%n)%n;
	gap+=h;
	floatbtmp=a.getview();
	//std::cout<<floatbtmp.num.size()<<"\n";
	ll diff=gap+a.ZERO*Blen;
	ll mindiff=(l*3+10-(int)floatbtmp.num.size())*Blen+diff%Blen;
	//std::cout<<diff<<" "<<mindiff<<"\n";
	int hiddenzero;
	if(diff>mindiff){floatbtmp.mul2pow(mindiff);hiddenzero=(diff-mindiff)/Blen;}
	else{floatbtmp.mul2pow(diff);hiddenzero=0;}
	//std::cout<<floatbtmp.num.size()<<"\n";
	a.NUM=floatbtmp.root(n,0,hiddenzero);
	ll aexp2=(a.EXP*Blen-gap)/n;
	int aexp=aexp2/Blen,v=aexp2%Blen;
	if(v<0){v+=Blen,aexp--;}
	if(v){a.NUM.mul2pow(v);}
	a.EXP=aexp;
	a.ZERO=0;
	a.LOW=0;
	a.setvalidlength(l);
	a.foldzero();
	return a;
}

integer floatb::stringtob28(const std::string&s,int&e1)
{
	integer tmp;
	int i=s.size()-1,end=std::max(-1,i-11);
	while(i>end&&s[i]!='e'){i--;}if(i==end){i=s.size();}
	int j=0;
	for(;j<i&&s[j]!='.';j++);if(j==i){j=-1;}
	e1=j>-1?-(i-j-1):0;
	if(i+1<s.size()){e1+=std::stoi(s.substr(i+1));}
	int k=i-8;
	for (; k> j&&k>0; k -= 8)
	{
		int a = s[k] - 48;
		for (int t = 1; t < 8; t++)
		{
			a *= 10;
			a += s[k + t] - 48;
		}
		tmp.num.push_back(a);
	}
	if(k>0)
	{
		k--;
		int a =s[k] - 48;
		int t=1;
		for (; t < j-k; t++)
		{
			a *= 10;
			a += s[k + t] - 48;
		}
		for (t++; t <9; t++)
		{
			a *= 10;
			a += s[k + t] - 48;
		}
		tmp.num.push_back(a);
		k-=8;
	}
	for (; k> 0; k -= 8)
	{
		int a = s[k] - 48;
		for (int t = 1; t < 8; t++)
		{
			a *= 10;
			a += s[k + t] - 48;
		}
		tmp.num.push_back(a);
	}
	k += 8;
	int a = 0, l = 0;
	if (s[0] == '-') { l = 1; tmp.sign = -1; }
	for (; l < k; l++)
	{
		if(l==j){continue;}
		a *= 10;
		a += s[l] - 48;
	}
	if (a || tmp.num.empty()) { tmp.num.push_back(a); }
	//for(auto&e:tmp.num){std::cout<<e<<" ";}std::cout<<"\n"<<exp*8<<"\n";
	tmp.e8tob28();
	return tmp;
}

floatb::floatb(const std::string&s)
{
	int e1;
	floatb tmp(stringtob28(s,e1),0);
	int l=std::max(minvalidlength,(int)tmp.NUM.num.size());
	tmp.setvalidlength(l);
	floatb h=integer(10);
	h.setvalidlength(l);
	*this=h.calpower(e1,l).multiply(tmp,l);
	foldzero();
}

void floatb::get()
{
	std::string s;
	std::cin>>s;
	*this=floatb(s);
}

std::string floatb::tostring(int m) const
{
	if(!NUM.num.back()){return std::string("0");}
    // 有一定误差
    integer::view aview = getview();
	integer tmp;
	tmp.num.reserve(aview.len+ZERO);
	tmp.num.assign(ZERO,0);
	for(int i=0;i<aview.len;i++){tmp.num.push_back(aview.ptr[i]);}
	aview=tmp;
    if (m > aview.len || m < 1) { m = aview.len; }
    int endexp = aview.len - m + EXP;
    int n = abs(endexp) * (log(Base) / log(e8));
    const bool flag = endexp < 0;
    if (flag) { n++; }
    else { n -= (n != 0); }
    int need = m + 2, b = 0;
    integer u = integer::shiftpow(e8, n, need, b);
    integer v;
    int vexp = EXP;
    if (flag)
    {
        v = integer::shiftmul(u, aview, need, b);
        vexp += b;
    }
    else
    {
        int gap = (int)u.num.size() - need;
        integer::view uview;
        if (gap > 0) { b += gap; uview = integer::view(u, gap); }
        else { uview = u; }
        gap = aview.len - need - uview.len;
        // 本来只要商 m 位的,商 need 保证 vexp<0
        vexp += gap;
        if (gap > 0) { floatbtmp = integer::view(aview, gap); }
        else
        {
            floatbtmp.num.reserve(-gap + aview.len);
            floatbtmp.num.assign(-gap, 0);
            for (int i = 0; i < aview.len; i++) { floatbtmp.num.push_back(aview.ptr[i]); }
        }
        v = integer::divide(floatbtmp, uview, floatbtmp);
        // 这里不能调用 floatb::divide, 完整算也还好
        vexp -= b;
    }
    if (vexp <= 0) { v.div2pow(-Blen * vexp); }
    else { std::cout << "print"; exit(0); }
    v.sign = NUM.sign;
    return v.tostring(m, flag ? -n : n);
}

void floatb::print(bool stdendl, int m) const
{
    std::cout << tostring(m);
    if (stdendl) { std::cout << '\n'; }
    else { std::cout << ' '; }
}
bool floatb::absbigger(const floatb&b)
{
	if(!NUM.num.back()){return 0;}
	int gap=validlength()+EXP-b.validlength()-b.EXP;
	if(gap){return gap>0;}
	integer::view  av=getview(),bv=b.getview();
	gap=bv.len-av.len;
	int i=av.len-1,end=std::max(-1,-gap-1);
	for(;i>end;i--)
	{
		if(av.ptr[i]!=bv.ptr[i+gap]){break;}
	}
	if(i==end){return av.len>bv.len;}
	return av.ptr[i]>bv.ptr[i+gap];
}
#endif
