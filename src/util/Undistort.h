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


#pragma once

#include <Eigen/Core>

#include "util/ImageAndExposure.h"
#include "util/NumType.h"


namespace dso
{

class Undistort {
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

    Undistort(Mat33f &K_, int wOrg_, int hOrg_, int w_, int h_);
    virtual ~Undistort();

    virtual void distortCoordinates(float* in_x, float* in_y, float* out_x,
                                    float* out_y, int n) const = 0;


    inline const Mat33f getK() const {
        return K;
    };
    inline const Eigen::Vector2i getSize() const {
        return Eigen::Vector2i(w, h);
    };
    inline const VecX getOriginalParameter() const {
        return parsOrg;
    };
    inline const Eigen::Vector2i getOriginalSize() const {
        return Eigen::Vector2i(wOrg, hOrg);
    };

    template<typename T> ImageAndExposure* undistort(ImageAndExposure* output) const;

    void loadPhotometricCalibration(std::string file, std::string noiseImage,
                                    std::string vignetteImage);

protected:
    int wOrg, hOrg, w, h, wUp, hUp;
    int upsampleUndistFactor;
    Mat33f K;
    VecX parsOrg;

    float* remapX;
    float* remapY;

    void applyBlurNoise(float* img) const;

    void makeOptimalK_crop();
    void makeOptimalK_full();

    int readFromFile(const char* configFileName, int nPars);
};


Undistort* getUndistorterFromFile(std::string configFilename);


class UndistortFOV : public Undistort
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

    UndistortFOV(Mat33f &K_, int wOrg_, int hOrg_, int w_, int h_,
                 float fx_, float fy_,
                 float cx_, float cy_, float dist_);
    ~UndistortFOV();
    void distortCoordinates(float* in_x, float* in_y, float* out_x, float* out_y, int n) const;

private:
    float fx;
    float fy;
    float cx;
    float cy;
    float dist;
};

class UndistortRadTan : public Undistort
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

    UndistortRadTan(Mat33f &K_, int wOrg_, int hOrg_, int w_, int h_,
                    float fx_, float fy_,
                    float cx_, float cy_,
                    float k1_, float k2_,
                    float r1_, float r2_);
    ~UndistortRadTan();
    void distortCoordinates(float* in_x, float* in_y, float* out_x, float* out_y,
                            int n) const;
private:
    float fx;
    float fy;
    float cx;
    float cy;
    float k1;
    float k2;
    float r1;
    float r2;
};
}

