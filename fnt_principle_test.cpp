/*
ai上找的strassen fft的核心费马数变换实现
自己推导的数学原理,可能不严谨或者有错,但是这个测试是对的
2^s作为limb很快,另外依赖mod x^n+1和 mod x^n-1实现crt递归，加权的负循环卷积实现规模折半
以下是推导
n=2^a,在mod=2^n+1下(不要求模质数)
2^n=-1,2^(2n)=1
选w0=2,长度2n,w=w0^t
证明dft公式存在即证明1+w1+w2...+w^(2n-1)=0
左边=(1+w^(2^a))*(1+w^(2^(a-1)))...*(1+w)
因为w^(2n)=1,设w阶2^m,1<=m<=a+1
w^(2m)=1得(w^(2^(m-1)-1)*(w^(2^(m-1)+1)=0得w^(2^(m-1))+1=0
故左边=0
另外w=2^t在2^n+1下逆元=2^(2n-t)=mod-2^(n-t)
在reverse需要2n在2^n+1下逆元
2逆元是(2^n+1)/2
而n=2^a,(2^a)*x+(2^(2^a)+1)*y=1,可以x=-2^(2^a-a),y=1,
所以n逆元是mod-2^(2^a-a)
*/
//请使用integer_old.h运行
#include"integer.h"
using namespace std;
integer N,MOD;
void ntt_test(vector<integer>&x,bool reverse)
{
    for (int i = 1, j = 0; i < x.size(); i++)
    {
        int bit = x.size() >> 1;
        for (; j & bit; j ^= bit, bit >>= 1); j ^= bit;
        if (i > j) { swap(x[i], x[j]); }
    }
    for(int len=2;len<=x.size();len<<=1)
    {
        integer w=1,w0=reverse?power(2,x.size()/len,MOD):MOD-power(2, N- x.size()/ len, MOD);
        for(int k=0;k<len/2;k++)
        {
            for(int i=0;i<x.size();i+=len)
            {
                integer  u=x[i+k],v=x[i+k+len/2]*w%MOD;
                x[i+k]=(u+v)%MOD,x[i+k+len/2]=(u+MOD-v)%MOD;
            }
            w=w*w0%MOD;
        }
    }
    if (reverse) 
    {   
        integer m=N;int  cnt=ctz(m);//ctz会改变传入的参数
        integer inv = (MOD+1)/2*(MOD-power(2,N-cnt,MOD))%MOD; 
        for (auto& e : x) { e = e * inv % MOD; } 
    }
}
int main()
{
   int n=11,m=1<<n;//n不要大于12,integer很慢的
   N=m/2,MOD=power(2,N)+1;
   vector<integer> x;
   x.assign(m,MOD);//预分配大小
   for(int i=0;i<m;i++){x[i]=integer(i)%MOD;}//不要破坏预分配的大小
   ntt_test(x,0);
   cout<<"MOD is ";MOD.print();
   ntt_test(x,1);
   n=1,m=0;
   for(auto&e:x){n&=(e.num[0]==(m++));}
   std::cout<<(n?"passed":"failed");
}
