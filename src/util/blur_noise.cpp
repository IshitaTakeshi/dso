#include "util/globalFuncs.h"
#include "util/blur_noise.h"

namespace dso {

BlurNoise::BlurNoise(const int w_, const int h_,
                     const int noise_grid_size_, const float noise_variance) :
                     w(w_), h(h_), noise_grid_size(noise_grid_size_) {
    assert(noise_variance >= 0);

    const int numnoise = (noise_grid_size+8) * (noise_grid_size+8);
    // Use smart pointer or
    for(int i=0; i<numnoise; i++) {
        noiseMapX.push_back(noise_variance * (rand()/(float)RAND_MAX));
        noiseMapY.push_back(noise_variance * (rand()/(float)RAND_MAX));
    }
}


void BlurNoise::apply(float* img) const {

    float blutTmp[w*h];

    float gaussMap[1000];
    for(int i=0; i<1000; i++) {
        gaussMap[i] = expf((float)(-i*i/(100.0*100.0)));
    }

    // x-blur.
    for(int y=0; y<h; y++) {
        for(int x=0; x<w; x++) {
            float xBlur = getInterpolatedElement11BiCub(
                noiseMapX,
                4+(x/(float)w)*noise_grid_size,
                4+(y/(float)h)*noise_grid_size,
                noise_grid_size+8
            );

            if(xBlur < 0.01) xBlur=0.01;


            int kernelSize = 1 + (int)(1.0f+xBlur*1.5);
            float sumW=0;
            float sumCW=0;
            for(int dx=0; dx <= kernelSize; dx++) {
                int gmid = 100.0f*dx/xBlur + 0.5f;
                if(gmid > 900) gmid = 900;
                float gw = gaussMap[gmid];

                if(x+dx>0 && x+dx<w) {
                    sumW += gw;
                    sumCW += gw * img[x+dx+y*this->w];
                }

                if(x-dx>0 && x-dx<w && dx!=0) {
                    sumW += gw;
                    sumCW += gw * img[x-dx+y*this->w];
                }
            }

            blutTmp[x+y*this->w] = sumCW / sumW;
        }
    }

    // y-blur.
    for(int x=0; x<w; x++) {
        for(int y=0; y<h; y++)
        {
            float yBlur = getInterpolatedElement11BiCub(
                              noiseMapY,
                              4+(x/(float)w)*noise_grid_size,
                              4+(y/(float)h)*noise_grid_size,
                              noise_grid_size+8);

            if(yBlur < 0.01) yBlur=0.01;

            int kernelSize = 1 + (int)(0.9f+yBlur*2.5);
            float sumW=0;
            float sumCW=0;
            for(int dy=0; dy <= kernelSize; dy++) {
                int gmid = 100.0f*dy/yBlur + 0.5f;
                if(gmid > 900 ) gmid = 900;
                float gw = gaussMap[gmid];

                if(y+dy>0 && y+dy<h) {
                    sumW += gw;
                    sumCW += gw * blutTmp[x+(y+dy)*this->w];
                }

                if(y-dy>0 && y-dy<h && dy!=0)
                {
                    sumW += gw;
                    sumCW += gw * blutTmp[x+(y-dy)*this->w];
                }
            }
            img[x+y*this->w] = sumCW / sumW;
        }
    }
}

}
