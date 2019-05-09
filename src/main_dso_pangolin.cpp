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



#include <thread>
#include <locale.h>
#include <signal.h>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>

#include "IOWrapper/Output3DWrapper.h"
#include "IOWrapper/ImageDisplay.h"

#include <boost/thread.hpp>
#include "util/settings.h"
#include "util/globalFuncs.h"
#include "util/DatasetReader.h"
#include "util/globalCalib.h"

#include "util/NumType.h"
#include "FullSystem/FullSystem.h"
#include "OptimizationBackend/MatrixAccumulators.h"
#include "FullSystem/PixelSelector2.h"

#include "IOWrapper/Pangolin/PangolinDSOViewer.h"
#include "IOWrapper/OutputWrapper/SampleOutputWrapper.h"


using namespace dso;


void my_exit_handler(int s)
{
    printf("Caught signal %d\n",s);
    exit(1);
}


void exitThread()
{
    struct sigaction sigIntHandler;
    sigIntHandler.sa_handler = my_exit_handler;
    sigemptyset(&sigIntHandler.sa_mask);
    sigIntHandler.sa_flags = 0;
    sigaction(SIGINT, &sigIntHandler, NULL);

    while(true) pause();
}


void parseArgument(char* arg,
                   std::string &vignette, std::string &gammaCalib,
                   std::string &source, std::string &calib,
                   float &playbackSpeed, bool &useSampleOutput) {
    int option;
    float foption;
    char buf[1000];

    if(1==sscanf(arg,"sampleoutput=%d",&option))
    {
        if(option==1)
        {
            useSampleOutput = true;
            printf("USING SAMPLE OUTPUT WRAPPER!\n");
        }
        return;
    }

    if(1==sscanf(arg,"quiet=%d",&option))
    {
        if(option==1)
        {
            setting_debugout_runquiet = true;
            printf("QUIET MODE, I'll shut up!\n");
        }
        return;
    }

    if(1==sscanf(arg,"rec=%d",&option))
    {
        if(option==0)
        {
            disableReconfigure = true;
            printf("DISABLE RECONFIGURE!\n");
        }
        return;
    }

    if(1==sscanf(arg,"noros=%d",&option))
    {
        if(option==1)
        {
            disableReconfigure = true;
            printf("DISABLE ROS (AND RECONFIGURE)!\n");
        }
        return;
    }

    if(1==sscanf(arg,"nomt=%d",&option))
    {
        if(option==1)
        {
            multiThreading = false;
            printf("NO MultiThreading!\n");
        }
        return;
    }

    if(1==sscanf(arg,"files=%s",buf))
    {
        source = buf;
        printf("loading data from %s!\n", source.c_str());
        return;
    }

    if(1==sscanf(arg,"calib=%s",buf))
    {
        calib = buf;
        printf("loading calibration from %s!\n", calib.c_str());
        return;
    }

    if(1==sscanf(arg,"vignette=%s",buf))
    {
        vignette = buf;
        printf("loading vignette from %s!\n", vignette.c_str());
        return;
    }

    if(1==sscanf(arg,"gamma=%s",buf))
    {
        gammaCalib = buf;
        printf("loading gammaCalib from %s!\n", gammaCalib.c_str());
        return;
    }

    if(1==sscanf(arg,"speed=%f",&foption))
    {
        playbackSpeed = foption;
        printf("PLAYBACK SPEED %f!\n", playbackSpeed);
        return;
    }

    if(1==sscanf(arg,"save=%d",&option))
    {
        if(option==1)
        {
            debugSaveImages = true;
            if(42==system("rm -rf images_out"))
                printf("system call returned 42 - what are the odds?. This is only here to shut up the compiler.\n");
            if(42==system("mkdir images_out"))
                printf("system call returned 42 - what are the odds?. This is only here to shut up the compiler.\n");
            if(42==system("rm -rf images_out"))
                printf("system call returned 42 - what are the odds?. This is only here to shut up the compiler.\n");
            if(42==system("mkdir images_out"))
                printf("system call returned 42 - what are the odds?. This is only here to shut up the compiler.\n");
            printf("SAVE IMAGES!\n");
        }
        return;
    }

    if(1==sscanf(arg,"mode=%d",&option))
    {

        if(option==0)
        {
            printf("PHOTOMETRIC MODE WITH CALIBRATION!\n");
        }
        if(option==1)
        {
            printf("PHOTOMETRIC MODE WITHOUT CALIBRATION!\n");
            setting_photometricCalibration = 0;
            setting_affineOptModeA = 0; //-1: fix. >=0: optimize (with prior, if > 0).
            setting_affineOptModeB = 0; //-1: fix. >=0: optimize (with prior, if > 0).
        }
        if(option==2)
        {
            printf("PHOTOMETRIC MODE WITH PERFECT IMAGES!\n");
            setting_photometricCalibration = 0;
            setting_affineOptModeA = -1; //-1: fix. >=0: optimize (with prior, if > 0).
            setting_affineOptModeB = -1; //-1: fix. >=0: optimize (with prior, if > 0).
            setting_minGradHistAdd=3;
        }
        return;
    }

    printf("could not parse argument \"%s\"!!!!\n", arg);
}


int main( int argc, char** argv )
{
    std::string vignette = "";
    std::string gammaCalib = "";
    std::string source = "";
    std::string calib = "";
    float playbackSpeed = 0;
    bool useSampleOutput = false;

    //setlocale(LC_ALL, "");
    for(int i=1; i<argc; i++) {
        parseArgument(argv[i],
                      vignette, gammaCalib, source, calib,
                      playbackSpeed, useSampleOutput);
    }

    std::cout << "vignette: " << vignette << std::endl;
    std::cout << "gammaCalib: " << gammaCalib << std::endl;
    std::cout << "source: " << source << std::endl;
    std::cout << "calib: " << calib << std::endl;

    // hook crtl+C.
    boost::thread exThread = boost::thread(exitThread);

    Undistort *undistort = getUndistorterFromFile(calib);
    if(undistort == 0) {
        exit(-1);
    }

    // TODO make photometricUndist independent from this class
    PhotometricUndistorter *photometricUndist = new PhotometricUndistorter(
        gammaCalib, "", vignette,
        undistort->getOriginalSize()[0], undistort->getOriginalSize()[1]);
    ImageFolderReader* reader = new ImageFolderReader(source, undistort, photometricUndist);
    reader->setGlobalCalibration();

    if(setting_photometricCalibration > 0 && photometricUndist->getG() == 0) {
        printf(
          "ERROR: dont't have photometric calibation. "
          "Need to use commandline options mode=1 or mode=2 "
        );
        exit(1);
    }

    FullSystem* fullSystem = new FullSystem(photometricUndist->getG(),
                                            undistort->getK(),
                                            playbackSpeed);

    IOWrap::PangolinDSOViewer* viewer = new IOWrap::PangolinDSOViewer(wG[0], hG[0], false);
    fullSystem->outputWrapper.push_back(viewer);

    if(useSampleOutput) {
        fullSystem->outputWrapper.push_back(new IOWrap::SampleOutputWrapper());
    }

    // to make MacOS happy: run this in dedicated thread -- and use this one to run the GUI.
    std::thread runthread([&]() {
        for(int i=0; i < reader->getNumImages(); i+=1) {
            ImageAndExposure* img = reader->getImage(i);

            fullSystem->addActiveFrame(img, i);

            delete img;

            if(fullSystem->initFailed || setting_fullResetRequested) {
                if(i < 250 || setting_fullResetRequested) {
                    printf("RESETTING!\n");

                    std::vector<IOWrap::Output3DWrapper*> wraps = fullSystem->outputWrapper;
                    delete fullSystem;

                    for(IOWrap::Output3DWrapper* ow : wraps) {
                        ow->reset();
                    }

                    fullSystem = new FullSystem(photometricUndist->getG(),
                                                undistort->getK(),
                                                playbackSpeed);
                    fullSystem->outputWrapper = wraps;

                    setting_fullResetRequested=false;
                }
            }

            if(fullSystem->isLost) {
                printf("LOST!!\n");
                break;
            }

        }
        fullSystem->blockUntilMappingIsFinished();

        fullSystem->printResult("result.txt");
    });

    viewer->run();

    runthread.join();

    for(IOWrap::Output3DWrapper* ow : fullSystem->outputWrapper) {
        ow->join();
        delete ow;
    }

    printf("DELETE FULLSYSTEM!\n");
    delete undistort;
    delete fullSystem;

    printf("DELETE READER!\n");
    delete reader;

    printf("EXIT NOW!\n");
    return 0;
}
