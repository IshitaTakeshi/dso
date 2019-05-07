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
#include "util/settings.h"
#include "util/globalFuncs.h"
#include "util/globalCalib.h"

#include <sstream>
#include <fstream>
#include <dirent.h>
#include <algorithm>

#include "util/Undistort.h"
#include "IOWrapper/ImageRW.h"

#include <boost/thread.hpp>

using namespace dso;



inline std::vector<std::string> getdir(std::string dirname) {
    std::vector<std::string> filenames;
    DIR *dp;
    struct dirent *dirp;

    if((dp = opendir(dirname.c_str())) == NULL) {
        return std::vector<std::string>();  // return empty list
    }

    while ((dirp = readdir(dp)) != NULL) {
        std::string name = std::string(dirp->d_name);

        if(name != "." && name != "..") {
            filenames.push_back(name);
        }
    }

    closedir(dp);

    for(auto &filename : filenames) {
        std::cout << "filename " << filename << std::endl;
    }

    std::sort(filenames.begin(), filenames.end());

    if(dirname[dirname.length() - 1] != '/') {
        dirname = dirname + "/";
    }
    for(std::string &filename : filenames) {
        filename = dirname + filename;
    }

    return filenames;
}


struct PrepImageItem
{
    int id;
    bool isQueud;
    ImageAndExposure* pt;

    inline PrepImageItem(int _id)
    {
        id=_id;
        isQueud = false;
        pt=0;
    }

    inline void release()
    {
        if(pt!=0) delete pt;
        pt=0;
    }
};


class ImageFolderReader
{
public:
    ImageFolderReader(std::string path, std::string calibFile,
                      std::string gammaFile, std::string vignetteFile) :
        undistort(Undistort::getUndistorterForFile(calibFile, gammaFile, vignetteFile)),
        filenames(getdir(path)) {

        if(filenames.size() == 0) {
            // TODO raise some exception
        }
        printf("ImageFolderReader: got %d files in %s!\n", (int)filenames.size(),
               path.c_str());

    }
    ~ImageFolderReader() {
        delete undistort;
    };

    Eigen::VectorXf getOriginalCalib() {
        return undistort->getOriginalParameter().cast<float>();
    }

    Eigen::Vector2i getOriginalDimensions() {
        return undistort->getOriginalSize();
    }

    void getCalibMono(Eigen::Matrix3f &K, int &w, int &h)
    {
        K = undistort->getK().cast<float>();
        w = undistort->getSize()[0];
        h = undistort->getSize()[1];
    }

    void setGlobalCalibration() {
        int w_out, h_out;
        Eigen::Matrix3f K;
        getCalibMono(K, w_out, h_out);
        setGlobalCalib(w_out, h_out, K);
    }

    int getNumImages() {
        return filenames.size();
    }

    ImageAndExposure* getImage(int id) {
        MinimalImageB* minimg = getImageRaw_internal(id);
        ImageAndExposure* ret2 = undistort->undistort<unsigned char>(minimg, 1.0f);
        delete minimg;
        return ret2;
    }

    inline float* getPhotometricGamma() {
        if(undistort==0 || undistort->photometricUndist==0) return 0;
        return undistort->photometricUndist->getG();
    }

private:
    std::vector<ImageAndExposure*> preloadedImages;
    const std::vector<std::string> filenames;
    // undistorter. [0] always exists, [1-2] only when MT is enabled.
    const Undistort* undistort;

    MinimalImageB* getImageRaw_internal(int id) {
        // CHANGE FOR ZIP FILE
        return IOWrap::readImageBW_8U(filenames[id]);
    }
};
