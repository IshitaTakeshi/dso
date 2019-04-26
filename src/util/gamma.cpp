#include "util/gamma.h"
#include <cassert>

namespace dso {

Gamma::Gamma(float *gammaInverse) : isAvailable(gammaInverse != 0) {
    if(!isAvailable) {
        return;
    }

    setGamma(gamma, gammaInverse);
}

void Gamma::setGamma(float* gamma, const float* gammaInv) {
    // invert.
    for(int i=1; i<255; i++) {
        // find val, such that Binv[val] = i.
        // I dont care about speed for this, so do it the stupid way.

        for(int s=1; s<255; s++) {
            if(gammaInv[s] <= i && i <= gammaInv[s+1]) {
                gamma[i] = s + (i - gammaInv[s]) / (gammaInv[s+1] - gammaInv[s]);
                break;
            }
        }
    }
    gamma[0] = 0;
    gamma[255] = 255;
}

float Gamma::get(int i) const {
    assert(isAvailable);
    return gamma[i];
}

}
