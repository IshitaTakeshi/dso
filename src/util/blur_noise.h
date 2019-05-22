#include <vector>

namespace dso {

class BlurNoise {
 public:
  BlurNoise(const int w_, const int h_,
            const int noise_grid_size_, const float noise_variance);
  void apply(float* img) const;
 private:
  std::vector<float> noiseMapX;
  std::vector<float> noiseMapY;
  const int w;
  const int h;
  const int noise_grid_size;
};

}
