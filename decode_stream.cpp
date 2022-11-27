#include <cmath>
#include <execution>
#include <fstream>
#include <iostream>
#include <random>
#include <unordered_map>
#include <sstream>
#include <string>

#include <dlib/optimization.h>

#include "asu.hpp"
#include "encode.hpp"
#include "geometry.hpp"
#include "partiality.hpp"
#include "wmath.hpp"
#include "wmath_optimisation.hpp"

using std::abs;
using std::array;
using std::cerr;
using std::cin;
using std::cout;
using std::endl;
using std::fill;
using std::fixed;
using std::get;
using std::getline;
using std::ifstream;
using std::isinf;
using std::isnan;
using std::isfinite;
using std::istream;
using std::lower_bound;
using std::make_tuple;
using std::max_element;
using std::nan;
using std::normal_distribution;
using std::numeric_limits;
using std::ofstream;
using std::random_device;
using std::round;
using std::unordered_map;
using std::unordered_set;
using std::setprecision;
using std::setw;
using std::sort;
using std::stod;
using std::streamsize;
using std::string;
using std::stringstream;
using std::swap;
using std::to_string;
using std::transform;
using std::tuple;
using std::vector;
using std::stoull;

using dlib::abs;
using dlib::cholesky_decomposition;
using dlib::identity_matrix;
using dlib::is_finite;
using dlib::length;
using dlib::length_squared;
using dlib::matrix;
using dlib::normalize;
using dlib::ones_matrix;
using dlib::squared;
using dlib::sum;
using dlib::tmp;
using dlib::trans;
using dlib::zeros_matrix;

using partiality::IDX;
using partiality::crystl;
using partiality::deserialize_crystls;
using partiality::deserialize_crystl;
using partiality::deserialize_sources;
using partiality::predict;
using partiality::predict_integrated;
using partiality::source;

using whash::patchmap;

using wmath::clip;
using wmath::count_stop_strategy;
using wmath::mean_variance;
using wmath::signum;
using wmath::pow;

using SYMMETRY::decode;
using SYMMETRY::get_point_group;
using SYMMETRY::reduce_encode;

constexpr double pi          = 3.14159265358979323846;

int main(int argc,char** argv) {
  for (size_t counter=0;cin;++counter) {
    auto sources = deserialize_sources(cin);
    if (!cin) break;
    for (auto it=sources.begin();it!=sources.end();++it) {
      cout << "> " << it->flx << " "
           << it->kin(0) << " "
           << it->kin(1) << " "
           << it->kin(2) << endl;
      cout << setw(16) << it->S12(0,0) << " "
           << setw(16) << it->S12(0,1) << " "
           << setw(16) << it->S12(0,2) << endl
           << setw(17)                 << " " 
           << setw(16) << it->S12(1,1) << " "
           << setw(16) << it->S12(1,2) << endl
           << setw(17)                 << " " 
           << setw(17)                 << " " 
           << setw(16) << it->S12(2,2) << endl;
    }
    auto crystl  = deserialize_crystl(cin);
    if (!cin) break;
    cout << "<" << endl;
    cout << setw(16) << crystl.R(0,0) << " "
         << setw(16) << crystl.R(0,1) << " "
         << setw(16) << crystl.R(0,2) << endl
         << setw(16) << crystl.R(1,0) << " "
         << setw(16) << crystl.R(1,1) << " "
         << setw(16) << crystl.R(1,2) << endl
         << setw(16) << crystl.R(2,0) << " "
         << setw(16) << crystl.R(2,1) << " "
         << setw(16) << crystl.R(2,2) << endl;
    cout << setw(16) << crystl.peak(0,0) << " "
         << setw(16) << crystl.peak(0,1) << " "
         << setw(16) << crystl.peak(0,2) << endl
         << setw(17)                     << " " 
         << setw(16) << crystl.peak(1,1) << " "
         << setw(16) << crystl.peak(1,2) << endl
         << setw(17)                     << " " 
         << setw(17)                     << " " 
         << setw(16) << crystl.peak(2,2) << endl;
    cout << crystl.mosaicity << " "
         << crystl.strain    << " "
         << crystl.a         << " "
         << crystl.b         << " ";
    double b = 1,c = 1;
    cin.read(reinterpret_cast<char*>(&b),sizeof(double));
    if (!cin) break;
    cout << b << " ";
    cin.read(reinterpret_cast<char*>(&c),sizeof(double));
    if (!cin) break;
    cout << c << endl;
    uint64_t m;
    cin.read(reinterpret_cast<char*>(&m),sizeof(uint64_t));
    int32_t h,k,l;
    float i,s;
    for (size_t j=0;j!=m;++j) {
      cin.read(reinterpret_cast<char*>(&h),sizeof(int32_t));
      cin.read(reinterpret_cast<char*>(&k),sizeof(int32_t));
      cin.read(reinterpret_cast<char*>(&l),sizeof(int32_t));
      cin.read(reinterpret_cast<char*>(&i),sizeof(float));
      cin.read(reinterpret_cast<char*>(&s),sizeof(float));
      //cout << h << " " << k << " " << l << " " << i << " " << s << endl;
      const IDX idx{h,k,l};
      const matrix<double,3,1> dhkl{1.0*h,1.0*k,1.0*l};
      const matrix<double,3,1> x = crystl.R*dhkl;
    }
  }
}
