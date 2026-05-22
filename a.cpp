
#include"integer.h"
#include<chrono>
using namespace std;
int main()
{
    integer n; n.get();
    // 只对 factor 计时
    auto start = std::chrono::high_resolution_clock::now();
    vector<integer>l;
    euler(n, &l);
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    cout << "\nfactor运行时间: " << duration.count() << " 微秒" << endl;
    for (auto& e : l) { e.print(); }
    return 0;
}//11111111111111111111111111111111111222233311123111111311
//1000000000000000000000000000156000000000000000000000000005643
//1111111111111111111111111111111111122223331112311111
//11111111111111111111111111111111111222233311123111111111111131
//1000000000000000000000000000156000000000000000000000011
