#include "vector.h"

using namespace NTL;
using namespace std;

namespace ntl {
/**
 * @brief Compute the norm of a vector
 */
double norm(Vec<ZZ> &v) {
    double res = 0;
    for (int i = 0; i < v.length(); i++) {
        res += conv<int>(v[i] * v[i]);
    }
    return sqrt(res);
}
} // namespace ntl
