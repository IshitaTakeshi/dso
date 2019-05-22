#include <cstring>


namespace dso {

class PhotometricUndistorter
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
    PhotometricUndistorter(std::string file, std::string noiseImage,
                           std::string vignetteImage, int w_, int h_);
    ~PhotometricUndistorter();

    // removes readout noise, and converts to irradiance.
    // affine normalizes values to 0 <= I < 256.
    // raw irradiance = a*I + b.
    template<typename T> float* processFrame(T* image_in) const;
    void unMapFloatImage(float* image);

    float* getG() {
        return G;
    };

private:
    float G[256*256];
    int GDepth;
    float* vignetteMap;
    float* vignetteMapInv;
    int w,h;
    bool valid;
};

}
