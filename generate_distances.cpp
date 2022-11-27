#include <algorithm>
#include <cmath>
#include <complex>
#include <cstddef>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <unordered_map>
#include <utility>

#include "wmath.hpp"
#include "patchmap.hpp"

using std::array;
using std::cerr;
using std::cin;
using std::complex;
using std::cout;
using std::endl;
using std::get;
using std::ifstream;
using std::isnan;
using std::normal_distribution;
using std::numeric_limits;
using std::ofstream;
using std::random_device;
using std::setw;
using std::sort;
using std::string;
using std::stringstream;
using std::tuple;
using std::uniform_real_distribution;
using std::unordered_map;
using std::vector;
using whash::patchmap;
using wmath::clip;
using wmath::mean_variance;
using wmath::pow;
using wmath::shl;
using wmath::shr;
using wmath::morton_encode;
using wmath::morton_decode;
using wmath::triangle_index_pack;
using wmath::triangle_index_unpack;

using namespace std::literals;

const double pi = 3.14159265358979323846;
const double  phi = (sqrt(5.0) + 1.0) * 0.5;
const double iphi = 1.0/phi;

// h k l i σ

struct dist_entry{
  double zscore = 0;
  double sumw  = 0;
  uint64_t n  = 0;
};

uint32_t pack(int32_t h,int32_t k, int32_t l){
  return
     shl(wmath::zigzag_encode(l),20u)
    ^shl(wmath::zigzag_encode(k),10u)
    ^shl(wmath::zigzag_encode(h),  0);
}

tuple<int32_t,int32_t,int32_t> unpack(uint32_t x){
  return {shr(x,20u)&(1023),shr(x,10)&(1023u),x&(1023u)};
}

typedef tuple<int16_t,int16_t,int16_t> millerindex;

int main(){
  patchmap<uint32_t,dist_entry> distances;
  vector<tuple<size_t,double,double>> data;
  size_t n = 0;
  for (string line;getline(cin,line);) { 
    stringstream ss(line);
    uint64_t n0;
    int32_t h,k,l;
    ss >> n0 >> h >> k >> l;
    if (!ss) {
      cerr << "could not read line:" << endl;
      cerr << line << endl;
      return 1;
    }
    if (ss.eof()) {
      data.clear();
      //cerr << '\r' << setw(16) << n << setw(16) << distances.size();
      ++n;
      continue;
    }
    double i0,v0,p0,e0,p,o,p1,e1;
    ss >> i0 >> v0 >> p0 >> e0 >> p >> o >> p1 >> e1;
    if (p0<=1e-8) continue;
    v0+=e0;
    i0/=p0;
    v0/=pow(p0,2);
    if (!ss) {
      cerr << "could not read line:" << endl;
      cerr << line << endl;
      return 1;
    }
    //const double s0 = p1/p0;
    //if (s0<0.25) continue;
    const double s0 = p0;
    for (size_t i=0;i!=data.size();++i) {
      const auto& [n1,i1,v1] = data[i];
      // x/s0 = y/s1 <-> x*s1/s0 = y <-> x = y*s0/s1
      //const double w = 1.0;
      //const double w = sqrt(pow(i0,2)+pow(i1,2));
      const double w = 1.0/(v0+v1);
      const uint32_t idx = triangle_index_pack(uint32_t(n1),uint32_t(n0));
      auto& entry = distances[idx];
      const double d = pow((i0-i1),2)/(v0+v1);
      entry.zscore += w*exp(2)*(1-exp(-d/exp(2)));
      entry.sumw   += w;
      entry.n      += 1;
    }
    data.emplace_back(n0,i0,v0);
  }
  //cerr << endl;
  for (auto it=distances.begin();it!=distances.end();++it) {
    const auto& [n0,n1] = triangle_index_unpack(it->first);
    const double& zscore = it->second.zscore;
    const double& sumw  = it->second.sumw;
    const uint32_t& n   = it->second.n;
    if (n<4) continue;
    cout << n0 << " " << n1 << " "
         << zscore/sumw << " "
         << n << " "
         << sumw << endl;
  }
}
// awk 'BEGIN{n=0}NF==4{delete i;delete v;n=0;next}{if ($7**2<1) next;n=$1;i[n]=$5/$7;v[n]=($6+$8)/$7**2;for (j in i) {if (j==n) continue; d[j][n]+=exp(2)*(1-exp(-(i[n]-i[j])**2/(v[n]+v[j])/exp(2)))/(v[n]+v[j]);s[j][n]+=1/(v[n]+v[j]);++num[j][n]};++n}END{for (j in d) for (n in d[j]) print j,n,d[j][n]/s[j][n],num[j][n]}'
