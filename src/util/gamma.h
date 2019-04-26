#ifndef UTIL_GAMMA_H
#define UTIL_GAMMA_H

namespace dso {

class Gamma {
 public:
    const bool isAvailable;
    Gamma(float *gammaInverse);
    float get(int i) const;

 private:
    float gamma[256];
    void setGamma(float* gamma, const float* gammaInv);
};

}

#endif /* UTIL_GAMMA_H */
