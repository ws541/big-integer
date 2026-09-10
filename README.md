# big-integer
# headonly  2k行 见integerb.h
# Base=2^28
# 支持快速乘法(karatsuba和FFT),高效除法(牛顿迭代法bunikel除法和原生除法),组合数学函数(阶乘、排列、组合计算),模运算系统(Montgomery模运算和通用模类),高精度开平方,高精度开n次方,最大公约数(hgcd、lehmer包含扩展欧几里得算法),模逆元,快速模幂运算,素数判定(Miller-Rabin和Lucas测试),jacobi,shanks模开方,质因数分解(Pollard Rho和二次筛法),欧拉函数,莫比乌斯函数计算以及原根求解,斐波那契数列,离散对数求解(bsgs+hellman+indexcalculus)

#  速度比专业库慢3到13倍,超大数字渐进行为良好，以下是部分测试结果
/*
测试代码见full_test.cpp
PS D:\cpp\.vscode\biginteger> ./a
=== integer Euclid vs GMP exGCD (Fibonacci) ===
n               size            integer(us)     GMP(us)         GMP快几倍
-----------------------------------------------------------------
n=100           3               49              71              0.69x
n=500           13              73              22              3.32x
n=2500          62              391             83              4.71x
n=12500         310             4241            799             5.31x
n=62500         1550            72910           14332           5.09x
n=312500                7749            490901          134024          3.66x

=== integer root7 vs GMP root7 (fac n) ===
n               size            integer(us)     GMP(us)         GMP快几倍
-----------------------------------------------------------------
n=100           19              117             24              4.88x
n=500           135             227             28              8.11x
n=2500          880             1214            154             7.88x
n=12500         5432            13645           1847            7.39x
n=62500         32342           81998           24718           3.32x

=== fac -> tostring -> construct (n=100000) ===
操作            integer(ms)     GMP(ms)         GMP快几倍
--------------------------------------------------------
fac             39              18              2.17x
tostring        210             66              3.18x
construct       65              31              2.10x
--------------------------------------------------------
total           314             115              2.73x
*/

#  数论方面函数因为小数字heavy慢,二次筛法和指标及算法进一步因为没有线性代数支持光滑界保守,差距可达几百倍
