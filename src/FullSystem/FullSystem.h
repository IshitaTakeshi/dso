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
#define MAX_ACTIVE_FRAMES 100

#include <deque>
#include <vector>

#include <iostream>
#include <fstream>

#include "util/globalCalib.h"
#include "util/NumType.h"
#include "util/gamma.h"
#include "util/FrameShell.h"
#include "util/IndexThreadReduce.h"
#include "OptimizationBackend/EnergyFunctional.h"
#include "FullSystem/HessianBlocks.h"
#include "FullSystem/PixelSelector2.h"
#include "FullSystem/Residuals.h"

namespace dso
{
namespace IOWrap
{
class Output3DWrapper;
}

class PixelSelector;
class PCSyntheticPoint;
class CoarseTracker;
struct FrameHessian;
struct PointHessian;
class CoarseInitializer;
struct ImmaturePointTemporaryResidual;
class CoarseDistanceMap;

class EnergyFunctional;

template<typename T> inline void deleteOut(std::vector<T*> &v, const int i)
{
    delete v[i];
    v[i] = v.back();
    v.pop_back();
}
template<typename T> inline void deleteOutPt(std::vector<T*> &v, const T* i)
{
    delete i;

    for(unsigned int k=0; k<v.size(); k++)
        if(v[k] == i)
        {
            v[k] = v.back();
            v.pop_back();
        }
}
template<typename T> inline void deleteOutOrder(std::vector<T*> &v,
        const int i)
{
    delete v[i];
    for(unsigned int k=i+1; k<v.size(); k++)
        v[k-1] = v[k];
    v.pop_back();
}
template<typename T> inline void deleteOutOrder(std::vector<T*> &v,
        const T* element)
{
    int i=-1;
    for(unsigned int k=0; k<v.size(); k++)
    {
        if(v[k] == element)
        {
            i=k;
            break;
        }
    }
    assert(i!=-1);

    for(unsigned int k=i+1; k<v.size(); k++)
        v[k-1] = v[k];
    v.pop_back();

    delete element;
}


inline bool eigenTestNan(const MatXX &m, std::string msg)
{
    bool foundNan = false;
    for(int y=0; y<m.rows(); y++)
        for(int x=0; x<m.cols(); x++)
        {
            if(!std::isfinite((double)m(y,x))) foundNan = true;
        }

    if(foundNan)
    {
        printf("NAN in %s:\n",msg.c_str());
        std::cout << m << "\n\n";
    }


    return foundNan;
}

void setPrecalcValues(std::vector<FrameHessian*> frameHessians, const Mat33f &K);

class FullSystem {
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    FullSystem(float *gammaInverse, const Mat33f &K);
    virtual ~FullSystem();

    // adds a new frame, and creates point & residual structs.
    void addActiveFrame(float* image, int id, const float exposure_time = 1.0);
    // marginalizes a frame. drops / marginalizes points & residuals.
    void marginalizeFrame(FrameHessian* frame);
    void blockUntilMappingIsFinished();

    float optimize(int mnumOptIts);

    void printResult(std::string file);


    // contains pointers to active frames

    std::vector<IOWrap::Output3DWrapper*> outputWrapper;

    bool isLost;
    bool initFailed;
    bool initialized;

private:
    const Gamma gamma;

    CameraParameters camera_parameters;

    // opt single point
    int optimizePoint(PointHessian* point, int minObs, bool flagOOB);
    PointHessian* optimizeImmaturePoint(ImmaturePoint* point, int minObs,
                                        ImmaturePointTemporaryResidual* residuals);

    double linAllPointSinle(PointHessian* point, float outlierTHSlack, bool plot);

    // mainPipelineFunctions
    Vec4 trackNewCoarse(FrameHessian* fh);
    void traceNewCoarse(FrameHessian* fh);
    void activatePoints();
    void activatePointsMT();
    void activatePointsOldFirst();
    void flagPointsForRemoval();
    void makeNewTraces(FrameHessian* newFrame, float* gtDepth);
    void initializeFromInitializer(FrameHessian* newFrame);
    void flagFramesForMarginalization();

    void removeOutliers();

    double linearizeAll(const std::vector<PointFrameResidual*> activeResiduals,
                        const Mat33f &K,
                        const bool fixLinearization);

    void linearizeAll_Reductor(const bool fixLinearization,
                               std::vector<PointFrameResidual*> toRemove,
                               const std::vector<PointFrameResidual*> activeResiduals,
                               const Mat33f &K,
                               Vec10* stats);
    void activatePointsMT_Reductor(std::vector<PointHessian*>* optimized,
                                   std::vector<ImmaturePoint*>* toOptimize,
                                   int min, int max, Vec10* stats, int tid);
    void printOptRes(const Vec3 &res, double resL, double resM, double resPrior,
                     double LExact, float a, float b);

    std::vector<VecX> getNullspaces(
        std::vector<VecX> &nullspaces_pose,
        std::vector<VecX> &nullspaces_scale);

    // =================== changed by tracker-thread. protected by trackMutex ============
    boost::mutex trackMutex;
    std::vector<FrameShell*> allFrameHistory;
    CoarseInitializer* coarseInitializer;
    Vec5 lastCoarseRMSE;

    // ================== changed by mapper-thread. protected by mapMutex ===============
    boost::mutex mapMutex;
    std::vector<FrameShell*> allKeyFramesHistory;

    EnergyFunctional* ef;
    IndexThreadReduce<Vec10> treadReduce;

    float* selectionMap;
    PixelSelector* pixelSelector;
    CoarseDistanceMap* coarseDistanceMap;

    // ONLY changed in marginalizeFrame and addFrame.
    std::vector<FrameHessian*> frameHessians;
    float currentMinActDist;

    // mutex etc. for tracker exchange.
    boost::mutex
    coarseTrackerSwapMutex;			// if tracker sees that there is a new reference, tracker locks [coarseTrackerSwapMutex] and swaps the two.
    CoarseTracker*
    coarseTracker_forNewKF;			// set as as reference. protected by [coarseTrackerSwapMutex].
    CoarseTracker*
    coarseTracker;					// always used to track new frames. protected by [trackMutex].
    float minIdJetVisTracker, maxIdJetVisTracker;
    float minIdJetVisDebug, maxIdJetVisDebug;

    VecC current_camera_parameters;



    // mutex for camToWorl's in shells (these are always in a good configuration).
    boost::mutex shellPoseMutex;



    /*
     * tracking always uses the newest KF as reference.
     *
     */

    void makeKeyFrame( FrameHessian* fh);
    void makeNonKeyFrame( FrameHessian* fh);
    void deliverTrackedFrame(FrameHessian* fh, bool needKF);
    void mappingLoop();

    // tracking / mapping synchronization. All protected by [trackMapSyncMutex].
    boost::mutex trackMapSyncMutex;
    boost::condition_variable trackedFrameSignal;
    boost::condition_variable mappedFrameSignal;
    std::deque<FrameHessian*> unmappedTrackedFrames;
    int needNewKFAfter;	// Otherwise, a new KF is *needed that has ID bigger than [needNewKFAfter]*.
    boost::thread mappingThread;
    bool runMapping;
    bool needToKetchupMapping;

    int lastRefStopID;
};
}

