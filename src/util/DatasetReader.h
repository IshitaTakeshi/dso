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

    std::sort(filenames.begin(), filenames.end());

    if(dirname[dirname.length() - 1] != '/') {
        dirname = dirname + "/";
    }
    for(std::string &filename : filenames) {
        filename = dirname + filename;
    }

    return filenames;
}


class ImageFolderReader {
public:
    //TODO remove undistort from class members
    ImageFolderReader(std::string path, Undistort* undistort_,
                      PhotometricUndistorter *photometricUndist_) :
        undistort(undistort_),
        photometricUndist(photometricUndist_),
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

    int getNumImages() {
        return filenames.size();
    }

    ImageAndExposure* getImage(int id) {
        MinimalImageB* minimg = getImageRaw_internal(id);

        if(minimg->w != undistort->getOriginalSize()[0] ||
           minimg->h != undistort->getOriginalSize()[1]) {
            printf("Undistort::undistort: wrong image size\n");
            exit(1);
        }

        ImageAndExposure* output = photometricUndist->processFrame<unsigned char>(
            minimg->data, 1.0f, 1.0f);
        ImageAndExposure* ret2 = undistort->undistort<unsigned char>(output);

        delete minimg;
        delete output;
        return ret2;
    }

private:
    const std::vector<std::string> filenames;
    // undistorter. [0] always exists, [1-2] only when MT is enabled.
    const Undistort* undistort;
    const PhotometricUndistorter* photometricUndist;

    MinimalImageB* getImageRaw_internal(int id) {
        // CHANGE FOR ZIP FILE
        return IOWrap::readImageBW_8U(filenames[id]);
    }
};
