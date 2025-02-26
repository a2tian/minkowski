#include "vector.h"

using namespace NTL;
using namespace std;

namespace ntl {
/**
 * @brief Compute the norm of a vector
 */
RR norm(Vec<ZZ> &v) {
    ZZ res;
    for (int i = 0; i < v.length(); i++) {
        res += v[i] * v[i];
    }
    return sqrt(conv<RR>(res));
}
}  // namespace ntl
