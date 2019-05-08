/**
* This file is part of DSO.
*
* Copyright 2016 Technical University of Munich and Intel.
* Developed by Jakob Engel <engelj at in dot tum dot de>,
* for more information see <http://vision.in.tum.de/dso>.
* If you use this code, please cite the respective publications as
* listed on the above website.
*
* DSO is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
*
* DSO is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with DSO. If not, see <http://www.gnu.org/licenses/>.
*/




#include <sstream>
#include <fstream>
#include <iostream>

#include <Eigen/Core>
#include <iterator>
#include "util/settings.h"
#include "util/globalFuncs.h"
#include "IOWrapper/ImageDisplay.h"
#include "IOWrapper/ImageRW.h"
#include "util/Undistort.h"
#include "util/camera_matrix.h"


namespace dso {


int readImageSize(int &w, int &h, const std::string &line) {
    if(std::sscanf(line.c_str(), "%d %d", &w, &h) == 2) {
        return 0;
    }
    return -1;
}


int isRelativeFormat(const float cx, const float cy) {
    return cx < 1 && cy < 1;
}


int isRadTanParameters(std::string line) {
    float ic[8];
    return std::sscanf(line.c_str(), "%f %f %f %f %f %f %f %f",
                       &ic[0], &ic[1], &ic[2], &ic[3],
                       &ic[4], &ic[5], &ic[6], &ic[7]) == 8;
}


int isFOVParameters(std::string line) {
    float ic[5];
    return std::sscanf(line.c_str(), "%f %f %f %f %f",
                       &ic[0], &ic[1], &ic[2], &ic[3], &ic[4]) == 5;
}


void relativeToAbsolute(float &fx, float &fy, float &cx, float &cy,
                        const int wOrg, const int hOrg) {
    // rescale and substract 0.5 offset.
    // the 0.5 is because I'm assuming the calibration is given such that the pixel at (0,0)
    // contains the integral over intensity over [0,0]-[1,1], whereas I assume the pixel (0,0)
    // to contain a sample of the intensity ot [0,0], which is best approximated by the integral over
    // [-0.5,-0.5]-[0.5,0.5]. Thus, the shift by -0.5.
    fx = fx * wOrg;
    fy = fy * hOrg;
    cx = cx * wOrg - 0.5;
    cy = cy * hOrg - 0.5;
}


// FIXME this is same as readFOVParameters_
int readOutputParameters_(float &fx, float &fy, float &cx, float &cy, float &r,
                          const std::string &line) {
    if(std::sscanf(line.c_str(), "%f %f %f %f %f",
                   &fx, &fy, &cx, &cy, &r) == 5) {
        return 0;
    }
    return -1;
}


int readOutputParameters(Mat33f &K, int w, int h, const std::string &line) {
    float fx, fy, cx, cy, r;

    std::cout << "called readOutputParameters" << std::endl;

    if(readOutputParameters_(fx, fy, cx, cy, r, line) != 0) {
        return -1;
    }

    if(!isRelativeFormat(cx, cy)) {
        return -1;
    }

    std::cout << "output " << std::endl;
    std::cout << "fx, fy, cx, cy, r" << " "
              << fx << " " << fy << " "
              << cx << " " << cy << " " << r << " " << std::endl;

    relativeToAbsolute(fx, fy, cx, cy, w, h);

    std::cout << "output " << std::endl;
    std::cout << "fx, fy, cx, cy, r" << " "
              << fx << " " << fy << " "
              << cx << " " << cy << " " << r << " " << std::endl;
    K = initializeCameraMatrix(fx, fy, cx, cy);
    return 0;
}


int readFOVParameters_(float &fx, float &fy, float &cx, float &cy, float &dist,
                       const std::string &line) {
    char buf[1000];
    snprintf(buf, 1000, "%%f %%f %%f %%f %%f");

    if(std::sscanf(line.c_str(), buf, &fx, &fy, &cx, &cy, &dist) == 5) {
        return 0;
    }
    return -1;
}


int readRadTanParameters_(float &fx, float &fy, float &cx, float &cy,
                          float &k1, float &k2, float &r1, float &r2,
                          const std::string &line) {
    char buf[1000];
    snprintf(buf, 1000, "%%f %%f %%f %%f %%f %%f %%f %%f %%f %%f");
    if(std::sscanf(line.c_str(), buf,
                   &fx, &fy, &cx, &cy, &k1, &k2, &r1, &r2) == 8) {
        return 0;
    }
    return -1;
}


int readFOVParameters(float &fx, float &fy, float &cx, float &cy, float &dist,
                      int wOrg, int hOrg,
                      const std::string &line) {
    if(readFOVParameters_(fx, fy, cx, cy, dist, line) != 0) {
        return -1;
    }
    if(isRelativeFormat(cx, cy)) {
        relativeToAbsolute(fx, fy, cx, cy, wOrg, hOrg);
    }
    return 0;
}


int readRadTanParameters(float &fx, float &fy, float &cx, float &cy,
                         float &k1, float &k2, float &r1, float &r2,
                         int wOrg, int hOrg,
                         const std::string &line) {
    if(readRadTanParameters_(fx, fy, cx, cy, k1, k2, r1, r2, line) != 0) {
        return -1;
    }
    if(isRelativeFormat(cx, cy)) {
        relativeToAbsolute(fx, fy, cx, cy, wOrg, hOrg);
    }
    return 0;
}


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


Undistort::~Undistort()
{
    if(remapX != 0) delete[] remapX;
    if(remapY != 0) delete[] remapY;
}


Undistort* getUndistorterForFile(std::string configFilename) {
    printf("Reading Calibration from file %s", configFilename.c_str());

    // read parameters
    std::ifstream infile(configFilename);
    assert(infile.good());

    std::string l1, l2, l3, l4;

    std::getline(infile, l1);
    std::getline(infile, l2);
    std::getline(infile, l3);
    std::getline(infile, l4);

    if (!infile.good()) {
        printf(" ... not found. Cannot operate without calibration, shutting down.\n");
        infile.close();
        return 0;
    }

    infile.close();

    Mat33f K;

    int wOrg, hOrg, w, h;

    // read image sizes first
    if(readImageSize(wOrg, hOrg, l2) != 0) {
        printf("Failed to read input image size\n");
        return 0;
    }

    if(readImageSize(w, h, l4) != 0) {
        printf("Failed to read output image size\n");
        return 0;
    }

    if(readOutputParameters(K, w, h, l3) != 0) {
        printf("Failed to read output camera parameters\n");
        return 0;
    }

    if(isRadTanParameters(l1)) {
        printf("found RadTan (OpenCV) camera model, building rectifier.\n");
        float fx, fy, cx, cy, k1, k2, r1, r2;
        readRadTanParameters(fx, fy, cx, cy, k1, k2, r1, r2, wOrg, hOrg, l1);
        return new UndistortRadTan(K, wOrg, hOrg, w, h, fx, fy, cx, cy, k1, k2, r1, r2);
    }

    if(isFOVParameters(l1)) {
        printf("found FOV camera model, building rectifier.\n");
        float fx, fy, cx, cy, dist;
        readFOVParameters(fx, fy, cx, cy, dist, wOrg, hOrg, l1);
        std::cout << "fx, fy, cx, cy, dist" << " "
                  << fx << " " << fy << " "
                  << cx << " " << cy << " " << dist << " " << std::endl;
        std::cout << "wOrg, hOrg, w, h" << " "
                  <<  wOrg << " " << hOrg << " "
                  << w << " " << h << " " << std::endl;
        std::cout << "K" << std::endl;
        std::cout << K << std::endl;
        return new UndistortFOV(K, wOrg, hOrg, w, h, fx, fy, cx, cy, dist);
    }

    printf("could not read calib file! exit.");
    return 0;
}



template<typename T>
ImageAndExposure* Undistort::undistort(ImageAndExposure* output) const {
    ImageAndExposure* result = new ImageAndExposure(w, h);
    output->copyMetaTo(*result);

    float* out_data = result->image;
    float* in_data = output->image;

    float* noiseMapX=0;
    float* noiseMapY=0;

    for(int idx = w*h-1; idx>=0; idx--) {
        // get interp. values
        float xx = remapX[idx];
        float yy = remapY[idx];

        if(xx<0)
            out_data[idx] = 0;
        else {
            // get integer and rational parts
            int xxi = xx;
            int yyi = yy;
            xx -= xxi;
            yy -= yyi;

            // get array base pointer
            const float* src = in_data + xxi + yyi * wOrg;

            // interpolate (bilinear)
            out_data[idx] = xx * yy * src[1+wOrg]
                          + (1-xx) * yy * src[wOrg]
                          + xx * (1-yy) * src[1]
                          + (1-xx) * (1-yy) * src[0];
        }
    }

    applyBlurNoise(result->image);

    return result;
}

template ImageAndExposure* Undistort::undistort<unsigned char> (ImageAndExposure* output) const;
template ImageAndExposure* Undistort::undistort<unsigned short> (ImageAndExposure* output) const;


void Undistort::applyBlurNoise(float* img) const
{
    if(benchmark_varBlurNoise==0) return;

    int numnoise=(benchmark_noiseGridsize+8)*(benchmark_noiseGridsize+8);
    float* noiseMapX = new float[numnoise];
    float* noiseMapY = new float[numnoise];
    float* blutTmp=new float[w*h];

    if(benchmark_varBlurNoise>0) {
        for(int i=0; i<numnoise; i++) {
            noiseMapX[i] =  benchmark_varBlurNoise  * (rand()/(float)RAND_MAX);
            noiseMapY[i] =  benchmark_varBlurNoise  * (rand()/(float)RAND_MAX);
        }
    }


    float gaussMap[1000];
    for(int i=0; i<1000; i++) {
        gaussMap[i] = expf((float)(-i*i/(100.0*100.0)));
    }

    // x-blur.
    for(int y=0; y<h; y++) {
        for(int x=0; x<w; x++) {
            float xBlur = getInterpolatedElement11BiCub(
                noiseMapX,
                4+(x/(float)w)*benchmark_noiseGridsize,
                4+(y/(float)h)*benchmark_noiseGridsize,
                benchmark_noiseGridsize+8
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
            float yBlur = getInterpolatedElement11BiCub(noiseMapY,
                          4+(x/(float)w)*benchmark_noiseGridsize,
                          4+(y/(float)h)*benchmark_noiseGridsize,
                          benchmark_noiseGridsize+8 );

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

    delete[] noiseMapX;
    delete[] noiseMapY;
}

// TODO take w and h as args or make them cost members
void makeRoundingResistant(float* remapX, float* remapY, int w, int h, int wOrg, int hOrg) {
    for(int y=0; y<h; y++) {
        for(int x=0; x<w; x++) {
            // make rounding resistant.
            float ix = remapX[y*w+x];
            float iy = remapY[y*w+x];

            if(ix == 0) ix = 0.001;
            if(iy == 0) iy = 0.001;
            if(ix == wOrg-1) ix = wOrg-1.001;
            if(iy == hOrg-1) ix = hOrg-1.001;

            if(ix > 0 && iy > 0 && ix < wOrg-1 && iy < hOrg-1) {
                remapX[y*w+x] = ix;
                remapY[y*w+x] = iy;
            } else {
                remapX[y*w+x] = -1;
                remapY[y*w+x] = -1;
            }
        }
    }
}


Undistort::Undistort(Mat33f &K_, int wOrg_, int hOrg_, int w_, int h_) :
    K(K_), wOrg(wOrg_), hOrg(hOrg_), w(w_), h(h_)
{
    remapX = new float[w*h];
    remapY = new float[w*h];
    for(int y=0; y<h; y++) {
        for(int x=0; x<w; x++) {
            remapX[y*w+x] = x;
            remapY[y*w+x] = y;
        }
    }
}


UndistortFOV::UndistortFOV(Mat33f &K_,
                           int wOrg_, int hOrg_, int w_, int h_,
                           float fx_, float fy_,
                           float cx_, float cy_, float dist_) :
    Undistort(K_, wOrg_, hOrg_, w_, h_),
    fx(fx_), fy(fy_), cx(cx_), cy(cy_), dist(dist_) {
    printf("Creating FOV undistorter\n");

    distortCoordinates(remapX, remapY, remapX, remapY, h*w);
    makeRoundingResistant(remapX, remapY, w, h, wOrg, hOrg);
}


UndistortFOV::~UndistortFOV() {}


void UndistortFOV::distortCoordinates(float* in_x, float* in_y, float* out_x,
                                      float* out_y, int n) const {
    float d2t = 2.0f * tan(dist / 2.0f);

    float ofx = K(0,0);
    float ofy = K(1,1);
    float ocx = K(0,2);
    float ocy = K(1,2);

    for(int i=0; i<n; i++) {
        float x = in_x[i];
        float y = in_y[i];
        float ix = (x - ocx) / ofx;
        float iy = (y - ocy) / ofy;

        float r = sqrtf(ix*ix + iy*iy);
        float fac = (r==0 || dist==0) ? 1 : atanf(r * d2t)/(dist*r);

        out_x[i] = fx*fac*ix+cx;
        out_y[i] = fy*fac*iy+cy;
    }
}

UndistortRadTan::UndistortRadTan(Mat33f &K_,
                                 int wOrg_, int hOrg_, int w_, int h_,
                                 float fx_, float fy_,
                                 float cx_, float cy_,
                                 float k1_, float k2_,
                                 float r1_, float r2_) :
    Undistort(K_, wOrg_, hOrg_, w_, h_),
    fx(fx_), fy(fy_), cx(cx_), cy(cy_), k1(k1_), k2(k2_), r1(r1_), r2(r2_) {

    printf("Creating RadTan undistorter\n");

    distortCoordinates(remapX, remapY, remapX, remapY, h*w);
    makeRoundingResistant(remapX, remapY, w, h, wOrg, hOrg);
}

UndistortRadTan::~UndistortRadTan() {}

void UndistortRadTan::distortCoordinates(float* in_x, float* in_y,
                                         float* out_x, float* out_y, int n) const {
    float ofx = K(0,0);
    float ofy = K(1,1);
    float ocx = K(0,2);
    float ocy = K(1,2);

    for(int i=0; i<n; i++) {
        float x = in_x[i];
        float y = in_y[i];

        // RADTAN
        float ix = (x - ocx) / ofx;
        float iy = (y - ocy) / ofy;
        float mx2_u = ix * ix;
        float my2_u = iy * iy;
        float mxy_u = ix * iy;
        float rho2_u = mx2_u+my2_u;
        float rad_dist_u = k1 * rho2_u + k2 * rho2_u * rho2_u;
        float x_dist = ix + ix * rad_dist_u + 2.0 * r1 * mxy_u + r2 *
                       (rho2_u + 2.0 * mx2_u);
        float y_dist = iy + iy * rad_dist_u + 2.0 * r2 * mxy_u + r1 *
                       (rho2_u + 2.0 * my2_u);
        out_x[i] = fx*x_dist+cx;
        out_y[i] = fy*y_dist+cy;
    }
}

}
