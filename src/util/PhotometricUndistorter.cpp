
#include <Eigen/Core>
#include "util/settings.h"
#include "util/PhotometricUndistorter.h"
#include "util/MinimalImage.h"
#include "util/ImageAndExposure.h"
#include "util/NumType.h"
#include "IOWrapper/ImageRW.h"

#include <iterator>
#include <fstream>


namespace dso {

PhotometricUndistorter::PhotometricUndistorter(
    std::string file,
    std::string noiseImage,
    std::string vignetteImage,
    int w_, int h_) {
    vignetteMap=0;
    vignetteMapInv=0;
    w = w_;
    h = h_;
    if(file=="" || vignetteImage=="") {
        printf("NO PHOTOMETRIC Calibration!\n");
    }

    // read G.
    std::ifstream f(file.c_str());
    printf("Reading Photometric Calibration from file %s\n",file.c_str());
    if (!f.good()) {
        printf("PhotometricUndistorter: Could not open file!\n");
        return;
    }

    {
        std::string line;
        std::getline(f, line);
        std::istringstream l1i(line);
        std::vector<float> Gvec = std::vector<float>(
            std::istream_iterator<float>
                                  (l1i), std::istream_iterator<float>() );

        GDepth = Gvec.size();

        if(GDepth < 256)
        {
            printf("PhotometricUndistorter: invalid format! "
                   "got %d entries in first line, expected at least 256!\n",
                   (int)Gvec.size());
            return;
        }


        for(int i=0; i<GDepth; i++) {
            G[i] = Gvec[i];
        }

        for(int i=0; i<GDepth-1; i++)
        {
            if(G[i+1] <= G[i])
            {
                printf("PhotometricUndistorter: G invalid! it has to be strictly increasing, but it isnt!\n");
                return;
            }
        }

        float min=G[0];
        float max=G[GDepth-1];
        for(int i=0; i<GDepth; i++) {
            G[i] = 255.0 * (G[i] - min) / (max-min);			// make it to 0..255 => 0..255.
        }
    }

    if(setting_photometricCalibration==0) {
        for(int i=0; i<GDepth; i++) G[i]=255.0f*i/(float)(GDepth-1);
    }

    printf("Reading Vignette Image from %s\n",vignetteImage.c_str());
    MinimalImage<unsigned short>* vm16 = IOWrap::readImageBW_16U(vignetteImage.c_str());
    MinimalImageB* vm8 = IOWrap::readImageBW_8U(vignetteImage.c_str());
    vignetteMap = new float[w*h];
    vignetteMapInv = new float[w*h];

    if(vm16 != 0)
    {
        if(vm16->w != w ||vm16->h != h)
        {
            printf("PhotometricUndistorter: Invalid vignette image size! got %d x %d, expected %d x %d\n",
                   vm16->w, vm16->h, w, h);
            if(vm16!=0) delete vm16;
            if(vm8!=0) delete vm8;
            return;
        }

        float maxV=0;
        for(int i=0; i<w*h; i++)
            if(vm16->at(i) > maxV) maxV = vm16->at(i);

        for(int i=0; i<w*h; i++)
            vignetteMap[i] = vm16->at(i) / maxV;
    }
    else if(vm8 != 0)
    {
        if(vm8->w != w ||vm8->h != h)
        {
            printf("PhotometricUndistorter: Invalid vignette image size! got %d x %d, expected %d x %d\n",
                   vm8->w, vm8->h, w, h);
            if(vm16!=0) delete vm16;
            if(vm8!=0) delete vm8;
            return;
        }

        float maxV=0;
        for(int i=0; i<w*h; i++)
            if(vm8->at(i) > maxV) maxV = vm8->at(i);

        for(int i=0; i<w*h; i++)
            vignetteMap[i] = vm8->at(i) / maxV;
    }
    else
    {
        printf("PhotometricUndistorter: Invalid vignette image\n");
        if(vm16!=0) delete vm16;
        if(vm8!=0) delete vm8;
        return;
    }

    if(vm16!=0) delete vm16;
    if(vm8!=0) delete vm8;


    for(int i=0; i<w*h; i++)
        vignetteMapInv[i] = 1.0f / vignetteMap[i];


    printf("Successfully read photometric calibration!\n");
}

PhotometricUndistorter::~PhotometricUndistorter() {
    if(vignetteMap != 0) delete[] vignetteMap;
    if(vignetteMapInv != 0) delete[] vignetteMapInv;
}


void PhotometricUndistorter::unMapFloatImage(float* image)
{
    int wh=w*h;
    for(int i=0; i<wh; i++)
    {
        float BinvC;
        float color = image[i];

        if(color < 1e-3)
            BinvC=0.0f;
        else if(color > GDepth-1.01f)
            BinvC=GDepth-1.1;
        else
        {
            int c = color;
            float a = color-c;
            BinvC=G[c]*(1-a) + G[c+1]*a;
        }

        float val = BinvC;
        if(val < 0) val = 0;
        image[i] = val;
    }
}

template<typename T> ImageAndExposure* PhotometricUndistorter::processFrame(
    T* image_in, float exposure_time, float factor) const {

    int wh=w*h;

    ImageAndExposure *output = new ImageAndExposure(w,h);
    float* data = output->image;
    assert(output->w == w && output->h == h);
    assert(data != 0);

    if(exposure_time <= 0 ||
       setting_photometricCalibration==0) { // disable full photometric calibration.
        for(int i=0; i<wh; i++) {
            data[i] = factor*image_in[i];
        }
    } else {
        for(int i=0; i<wh; i++) {
            data[i] = G[image_in[i]];
        }

        if(setting_photometricCalibration==2) {
            for(int i=0; i<wh; i++)
                data[i] *= vignetteMapInv[i];
        }
    }

    output->exposure_time = exposure_time;
    return output;
}

template ImageAndExposure* PhotometricUndistorter::processFrame<unsigned char>
(unsigned char* image_in, float exposure_time, float factor) const;
template ImageAndExposure* PhotometricUndistorter::processFrame<unsigned short>
(unsigned short* image_in, float exposure_time, float factor) const;

}
