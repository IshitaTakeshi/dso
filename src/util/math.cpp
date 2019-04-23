#include <algorithm>
#include <util/math.h>


using namespace std;

float clamp(float x, float bottom, float top) {
    return max(bottom, min(x, top));
}
