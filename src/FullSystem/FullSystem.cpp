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


/*
 * KFBuffer.cpp
 *
 *  Created on: Jan 7, 2014
 *      Author: engelj
 */

#include "stdio.h"
#include <Eigen/LU>
#include <algorithm>
#include <Eigen/SVD>
#include <Eigen/Eigenvalues>

#include "FullSystem/FullSystem.h"
#include "IOWrapper/ImageDisplay.h"
#include "util/globalFuncs.h"
#include "util/globalCalib.h"
#include "FullSystem/PixelSelector.h"
#include "FullSystem/PixelSelector2.h"
#include "FullSystem/ImmaturePoint.h"

#include "FullSystem/CoarseTracker.h"
#include "FullSystem/CoarseInitializer.h"

#include "OptimizationBackend/EnergyFunctional.h"
#include "OptimizationBackend/EnergyFunctionalStructs.h"

#include "IOWrapper/Output3DWrapper.h"

#include "util/camera_matrix.h"
#include "util/gamma.h"

#include <cmath>
#include <iostream>

namespace dso {


FullSystem::FullSystem(float *gammaInverse, const Mat33f &K, const float playbackSpeed) :
    gamma(Gamma(gammaInverse)), camera_parameters(CameraParameters(K)),
    linearizeOperation(playbackSpeed==0) {

    selectionMap = new float[wG[0]*hG[0]];

    coarseDistanceMap = new CoarseDistanceMap(wG[0], hG[0]);
    coarseTracker = new CoarseTracker(wG[0], hG[0]);
    coarseTracker_forNewKF = new CoarseTracker(wG[0], hG[0]);
    coarseInitializer = 0;
    pixelSelector = new PixelSelector(wG[0], hG[0]);

    current_camera_parameters = inv_scale_camera_parameters(camera_parameters.get());
    lastCoarseRMSE.setConstant(100);

    currentMinActDist=2;
    initialized=false;

    ef = new EnergyFunctional(K);
    ef->red = &this->treadReduce;

    isLost=false;
    initFailed=false;

    needNewKFAfter = -1;

    runMapping = true;
    mappingThread = boost::thread(&FullSystem::mappingLoop, this);
    lastRefStopID = 0;

    minIdJetVisDebug = -1;
    maxIdJetVisDebug = -1;
    minIdJetVisTracker = -1;
    maxIdJetVisTracker = -1;
}

FullSystem::~FullSystem() {
    blockUntilMappingIsFinished();

    delete[] selectionMap;

    for(FrameShell* s : allFrameHistory)
        delete s;
    for(FrameHessian* fh : unmappedTrackedFrames)
        delete fh;

    delete coarseDistanceMap;
    delete coarseTracker;
    delete coarseTracker_forNewKF;
    delete coarseInitializer;
    delete pixelSelector;
    delete ef;
}


void FullSystem::printResult(std::string file) {
    boost::unique_lock<boost::mutex> lock(trackMutex);
    boost::unique_lock<boost::mutex> crlock(shellPoseMutex);

    std::ofstream myfile;
    myfile.open (file.c_str());
    myfile << std::setprecision(15);
    myfile.close();
}


Vec4 FullSystem::trackNewCoarse(FrameHessian* fh) {
    assert(allFrameHistory.size() > 0);
    // set pose initialization.

    for(IOWrap::Output3DWrapper* ow : outputWrapper) {
        ow->pushLiveFrame(fh);
    }

    FrameHessian* lastF = coarseTracker->lastRef;

    AffLight aff_last_2_l = AffLight(0,0);

    std::vector<SE3,Eigen::aligned_allocator<SE3>> lastF_2_fh_tries;
    if(allFrameHistory.size() == 2) {
        for(unsigned int i=0; i<lastF_2_fh_tries.size(); i++) {
            lastF_2_fh_tries.push_back(SE3());
        }
    } else {
        FrameShell* slast = allFrameHistory[allFrameHistory.size()-2];
        FrameShell* sprelast = allFrameHistory[allFrameHistory.size()-3];
        SE3 slast_2_sprelast;
        SE3 lastF_2_slast;
        {   // lock on global pose consistency!
            boost::unique_lock<boost::mutex> crlock(shellPoseMutex);
            slast_2_sprelast = sprelast->camToWorld.inverse() * slast->camToWorld;
            lastF_2_slast = slast->camToWorld.inverse() * lastF->shell->camToWorld;
            aff_last_2_l = slast->aff_g2l;
        }
        SE3 fh_2_slast = slast_2_sprelast;// assumed to be the same as fh_2_slast.


        // get last delta-movement.
        // assume constant motion.
        lastF_2_fh_tries.push_back(
            fh_2_slast.inverse() * lastF_2_slast);
        // assume double motion (frame skipped)
        lastF_2_fh_tries.push_back(
            fh_2_slast.inverse() * fh_2_slast.inverse() * lastF_2_slast);
        lastF_2_fh_tries.push_back(SE3::exp(fh_2_slast.log()*0.5).inverse() *
                                   lastF_2_slast); // assume half motion.
        lastF_2_fh_tries.push_back(lastF_2_slast); // assume zero motion.
        lastF_2_fh_tries.push_back(SE3()); // assume zero motion FROM KF.

        // just try a TON of different initializations (all rotations). In the end,
        // if they don't work they will only be tried on the coarsest level, which is super fast anyway.
        // also, if tracking rails here we loose, so we really, really want to avoid that.
        for(float rotDelta=0.02; rotDelta < 0.05; rotDelta++)
        {
            // assume constant motion.
            lastF_2_fh_tries.push_back(
                fh_2_slast.inverse() * lastF_2_slast *
                    SE3(Sophus::Quaterniond(1,rotDelta,0,0), Vec3(0,0,0))
            );
            lastF_2_fh_tries.push_back(fh_2_slast.inverse() * lastF_2_slast * SE3(
                                           Sophus::Quaterniond(1,0,rotDelta,0), Vec3(0,0,
                                                   0)));			// assume constant motion.
            lastF_2_fh_tries.push_back(fh_2_slast.inverse() * lastF_2_slast * SE3(
                                           Sophus::Quaterniond(1,0,0,rotDelta), Vec3(0,0,
                                                   0)));			// assume constant motion.
            lastF_2_fh_tries.push_back(fh_2_slast.inverse() * lastF_2_slast * SE3(
                                           Sophus::Quaterniond(1,-rotDelta,0,0), Vec3(0,0,
                                                   0)));			// assume constant motion.
            lastF_2_fh_tries.push_back(fh_2_slast.inverse() * lastF_2_slast * SE3(
                                           Sophus::Quaterniond(1,0,-rotDelta,0), Vec3(0,0,
                                                   0)));			// assume constant motion.
            lastF_2_fh_tries.push_back(fh_2_slast.inverse() * lastF_2_slast * SE3(
                                           Sophus::Quaterniond(1,0,0,-rotDelta), Vec3(0,0,
                                                   0)));			// assume constant motion.
            lastF_2_fh_tries.push_back(fh_2_slast.inverse() * lastF_2_slast * SE3(
                                           Sophus::Quaterniond(1,rotDelta,rotDelta,0), Vec3(0,0,
                                                   0)));	// assume constant motion.
            lastF_2_fh_tries.push_back(fh_2_slast.inverse() * lastF_2_slast * SE3(
                                           Sophus::Quaterniond(1,0,rotDelta,rotDelta), Vec3(0,0,
                                                   0)));	// assume constant motion.
            lastF_2_fh_tries.push_back(fh_2_slast.inverse() * lastF_2_slast * SE3(
                                           Sophus::Quaterniond(1,rotDelta,0,rotDelta), Vec3(0,0,
                                                   0)));	// assume constant motion.
            lastF_2_fh_tries.push_back(fh_2_slast.inverse() * lastF_2_slast * SE3(
                                           Sophus::Quaterniond(1,-rotDelta,rotDelta,0), Vec3(0,0,
                                                   0)));	// assume constant motion.
            lastF_2_fh_tries.push_back(fh_2_slast.inverse() * lastF_2_slast * SE3(
                                           Sophus::Quaterniond(1,0,-rotDelta,rotDelta), Vec3(0,0,
                                                   0)));	// assume constant motion.
            lastF_2_fh_tries.push_back(fh_2_slast.inverse() * lastF_2_slast * SE3(
                                           Sophus::Quaterniond(1,-rotDelta,0,rotDelta), Vec3(0,0,
                                                   0)));	// assume constant motion.
            lastF_2_fh_tries.push_back(fh_2_slast.inverse() * lastF_2_slast * SE3(
                                           Sophus::Quaterniond(1,rotDelta,-rotDelta,0), Vec3(0,0,
                                                   0)));	// assume constant motion.
            lastF_2_fh_tries.push_back(fh_2_slast.inverse() * lastF_2_slast * SE3(
                                           Sophus::Quaterniond(1,0,rotDelta,-rotDelta), Vec3(0,0,
                                                   0)));	// assume constant motion.
            lastF_2_fh_tries.push_back(fh_2_slast.inverse() * lastF_2_slast * SE3(
                                           Sophus::Quaterniond(1,rotDelta,0,-rotDelta), Vec3(0,0,
                                                   0)));	// assume constant motion.
            lastF_2_fh_tries.push_back(fh_2_slast.inverse() * lastF_2_slast * SE3(
                                           Sophus::Quaterniond(1,-rotDelta,-rotDelta,0), Vec3(0,0,
                                                   0)));	// assume constant motion.
            lastF_2_fh_tries.push_back(fh_2_slast.inverse() * lastF_2_slast * SE3(
                                           Sophus::Quaterniond(1,0,-rotDelta,-rotDelta), Vec3(0,0,
                                                   0)));	// assume constant motion.
            lastF_2_fh_tries.push_back(fh_2_slast.inverse() * lastF_2_slast * SE3(
                                           Sophus::Quaterniond(1,-rotDelta,0,-rotDelta), Vec3(0,0,
                                                   0)));	// assume constant motion.
            lastF_2_fh_tries.push_back(fh_2_slast.inverse() * lastF_2_slast * SE3(
                                           Sophus::Quaterniond(1,-rotDelta,-rotDelta,-rotDelta), Vec3(0,0,
                                                   0)));	// assume constant motion.
            lastF_2_fh_tries.push_back(fh_2_slast.inverse() * lastF_2_slast * SE3(
                                           Sophus::Quaterniond(1,-rotDelta,-rotDelta,rotDelta), Vec3(0,0,
                                                   0)));	// assume constant motion.
            lastF_2_fh_tries.push_back(fh_2_slast.inverse() * lastF_2_slast * SE3(
                                           Sophus::Quaterniond(1,-rotDelta,rotDelta,-rotDelta), Vec3(0,0,
                                                   0)));	// assume constant motion.
            lastF_2_fh_tries.push_back(fh_2_slast.inverse() * lastF_2_slast * SE3(
                                           Sophus::Quaterniond(1,-rotDelta,rotDelta,rotDelta), Vec3(0,0,
                                                   0)));	// assume constant motion.
            lastF_2_fh_tries.push_back(fh_2_slast.inverse() * lastF_2_slast * SE3(
                                           Sophus::Quaterniond(1,rotDelta,-rotDelta,-rotDelta), Vec3(0,0,
                                                   0)));	// assume constant motion.
            lastF_2_fh_tries.push_back(fh_2_slast.inverse() * lastF_2_slast * SE3(
                                           Sophus::Quaterniond(1,rotDelta,-rotDelta,rotDelta), Vec3(0,0,
                                                   0)));	// assume constant motion.
            lastF_2_fh_tries.push_back(fh_2_slast.inverse() * lastF_2_slast * SE3(
                                           Sophus::Quaterniond(1,rotDelta,rotDelta,-rotDelta), Vec3(0,0,
                                                   0)));	// assume constant motion.
            lastF_2_fh_tries.push_back(fh_2_slast.inverse() * lastF_2_slast * SE3(
                                           Sophus::Quaterniond(1,rotDelta,rotDelta,rotDelta), Vec3(0,0,
                                                   0)));	// assume constant motion.
        }

        if(!slast->poseValid || !sprelast->poseValid || !lastF->shell->poseValid)
        {
            lastF_2_fh_tries.clear();
            lastF_2_fh_tries.push_back(SE3());
        }
    }


    Vec3 flowVecs = Vec3(100,100,100);
    SE3 lastF_2_fh = SE3();
    AffLight aff_g2l = AffLight(0,0);


    // as long as maxResForImmediateAccept is not reached, I'll continue through the options.
    // I'll keep track of the so-far best achieved residual for each level in achievedRes.
    // If on a coarse level, tracking is WORSE than achievedRes, we will not continue to save time.


    Vec5 achievedRes = Vec5::Constant(NAN);
    bool haveOneGood = false;
    int tryIterations=0;
    for(unsigned int i=0; i<lastF_2_fh_tries.size(); i++)
    {
        AffLight aff_g2l_this = aff_last_2_l;
        SE3 lastF_2_fh_this = lastF_2_fh_tries[i];
        bool trackingIsGood = coarseTracker->trackNewestCoarse(
                                  fh, lastF_2_fh_this, aff_g2l_this,
                                  pyrLevelsUsed-1,
                                  achievedRes);	// in each level has to be at least as good as the last try.
        tryIterations++;

        if(i != 0)
        {
            printf("RE-TRACK ATTEMPT %d with initOption %d and start-lvl %d (ab %f %f): %f %f %f %f %f -> %f %f %f %f %f \n",
                   i,
                   i, pyrLevelsUsed-1,
                   aff_g2l_this.a,aff_g2l_this.b,
                   achievedRes[0],
                   achievedRes[1],
                   achievedRes[2],
                   achievedRes[3],
                   achievedRes[4],
                   coarseTracker->lastResiduals[0],
                   coarseTracker->lastResiduals[1],
                   coarseTracker->lastResiduals[2],
                   coarseTracker->lastResiduals[3],
                   coarseTracker->lastResiduals[4]);
        }


        // do we have a new winner?
        if(trackingIsGood && std::isfinite((float)coarseTracker->lastResiduals[0])
                && !(coarseTracker->lastResiduals[0] >=  achievedRes[0]))
        {
            //printf("take over. minRes %f -> %f!\n", achievedRes[0], coarseTracker->lastResiduals[0]);
            flowVecs = coarseTracker->lastFlowIndicators;
            aff_g2l = aff_g2l_this;
            lastF_2_fh = lastF_2_fh_this;
            haveOneGood = true;
        }

        // take over achieved res (always).
        if(haveOneGood)
        {
            for(int i=0; i<5; i++)
            {
                if(!std::isfinite((float)achievedRes[i])
                        || achievedRes[i] >
                        coarseTracker->lastResiduals[i])	// take over if achievedRes is either bigger or NAN.
                    achievedRes[i] = coarseTracker->lastResiduals[i];
            }
        }


        if(haveOneGood &&  achievedRes[0] < lastCoarseRMSE[0]*setting_reTrackThreshold)
            break;

    }

    if(!haveOneGood)
    {
        printf("BIG ERROR! tracking failed entirely. Take predictred pose and hope we may somehow recover.\n");
        flowVecs = Vec3(0,0,0);
        aff_g2l = aff_last_2_l;
        lastF_2_fh = lastF_2_fh_tries[0];
    }

    lastCoarseRMSE = achievedRes;

    // no lock required, as fh is not used anywhere yet.
    fh->shell->camToTrackingRef = lastF_2_fh.inverse();
    fh->shell->trackingRef = lastF->shell;
    fh->shell->aff_g2l = aff_g2l;
    fh->shell->camToWorld = fh->shell->trackingRef->camToWorld *
                            fh->shell->camToTrackingRef;


    if(coarseTracker->firstCoarseRMSE < 0)
        coarseTracker->firstCoarseRMSE = achievedRes[0];

    if(!setting_debugout_runquiet)
        printf("Coarse Tracker tracked ab = %f %f (exp %f). Res %f!\n", aff_g2l.a,
               aff_g2l.b, fh->ab_exposure, achievedRes[0]);


    return Vec4(achievedRes[0], flowVecs[0], flowVecs[1], flowVecs[2]);
}

void FullSystem::traceNewCoarse(FrameHessian* fh)
{
    boost::unique_lock<boost::mutex> lock(mapMutex);
    Mat33f K = createCameraMatrixFromCalibHessian(camera_parameters);

    for(FrameHessian* host : frameHessians) {
        // go through all active frames
        SE3 hostToNew = fh->PRE_worldToCam * host->PRE_camToWorld;
        Mat33f KRKi = K * hostToNew.rotationMatrix().cast<float>() * K.inverse();
        Vec3f Kt = K * hostToNew.translation().cast<float>();

        Vec2f aff = AffLight::fromToVecExposure(
            host->ab_exposure, fh->ab_exposure,
            host->aff_g2l(), fh->aff_g2l()).cast<float>();

        for(ImmaturePoint* ph : host->immaturePoints) {
            ph->traceOn(KRKi, Kt, aff, fh->dI);
        }
    }
}


void FullSystem::activatePointsMT_Reductor(
    std::vector<PointHessian*>* optimized,
    std::vector<ImmaturePoint*>* toOptimize,
    int min, int max, Vec10* stats, int tid) {
    ImmaturePointTemporaryResidual* tr = new ImmaturePointTemporaryResidual[frameHessians.size()];
    for(int k=min; k<max; k++) {
        (*optimized)[k] = optimizeImmaturePoint((*toOptimize)[k], 1, tr);
    }
    delete[] tr;
}


void FullSystem::activatePointsMT()
{

    if(ef->nPoints < setting_desiredPointDensity*0.66)
        currentMinActDist -= 0.8;
    if(ef->nPoints < setting_desiredPointDensity*0.8)
        currentMinActDist -= 0.5;
    else if(ef->nPoints < setting_desiredPointDensity*0.9)
        currentMinActDist -= 0.2;
    else if(ef->nPoints < setting_desiredPointDensity)
        currentMinActDist -= 0.1;

    if(ef->nPoints > setting_desiredPointDensity*1.5)
        currentMinActDist += 0.8;
    if(ef->nPoints > setting_desiredPointDensity*1.3)
        currentMinActDist += 0.5;
    if(ef->nPoints > setting_desiredPointDensity*1.15)
        currentMinActDist += 0.2;
    if(ef->nPoints > setting_desiredPointDensity)
        currentMinActDist += 0.1;

    if(currentMinActDist < 0) currentMinActDist = 0;
    if(currentMinActDist > 4) currentMinActDist = 4;

    if(!setting_debugout_runquiet)
        printf("SPARSITY:  MinActDist %f (need %d points, have %d points)!\n",
               currentMinActDist, (int)(setting_desiredPointDensity), ef->nPoints);

    FrameHessian* newestHs = frameHessians.back();

    // make dist map.
    coarseDistanceMap->makeK(&camera_parameters);
    coarseDistanceMap->makeDistanceMap(frameHessians, newestHs);

    //coarseTracker->debugPlotDistMap("distMap");

    std::vector<ImmaturePoint*> toOptimize;
    toOptimize.reserve(20000);

    for(FrameHessian* host : frameHessians)		// go through all active frames
    {
        if(host == newestHs) continue;

        SE3 fhToNew = newestHs->PRE_worldToCam * host->PRE_camToWorld;
        Mat33f KRKi = (
            coarseDistanceMap->K[1]
            * fhToNew.rotationMatrix().cast<float>()
            * coarseDistanceMap->Ki[0]);
        Vec3f Kt = (coarseDistanceMap->K[1] * fhToNew.translation().cast<float>());

        for(unsigned int i=0; i<host->immaturePoints.size(); i+=1)
        {
            ImmaturePoint* ph = host->immaturePoints[i];
            ph->idxInImmaturePoints = i;

            // delete points that have never been traced successfully, or that are outlier on the last trace.
            if(!std::isfinite(ph->idepth_max) || ph->lastTraceStatus == IPS_OUTLIER)
            {
//				immature_invalid_deleted++;
                // remove point.
                delete ph;
                host->immaturePoints[i]=0;
                continue;
            }

            // can activate only if this is true.
            bool canActivate = (
                    ph->lastTraceStatus == IPS_GOOD ||
                    ph->lastTraceStatus == IPS_SKIPPED ||
                    ph->lastTraceStatus == IPS_BADCONDITION ||
                    ph->lastTraceStatus == IPS_OOB
                ) &&
                ph->lastTracePixelInterval < 8 &&
                ph->quality > setting_minTraceQuality &&
                (ph->idepth_max+ph->idepth_min) > 0;


            // if I cannot activate the point, skip it. Maybe also delete it.
            if(!canActivate)
            {
                // if point will be out afterwards, delete it instead.
                if(ph->host->flaggedForMarginalization || ph->lastTraceStatus == IPS_OOB)
                {
                    delete ph;
                    host->immaturePoints[i]=0;
                }
                continue;
            }

            // see if we need to activate point due to distance map.
            Vec3f ptp = KRKi * Vec3f(ph->u, ph->v, 1)
                      + Kt*(0.5f*(ph->idepth_max+ph->idepth_min));
            int u = ptp[0] / ptp[2] + 0.5f;
            int v = ptp[1] / ptp[2] + 0.5f;

            if((u > 0 && v > 0 && u < wG[1] && v < hG[1])) {
                float dist = coarseDistanceMap->fwdWarpedIDDistFinal[u+wG[1]*v] +
                             (ptp[0]-floorf((float)(ptp[0])));

                if(dist>=currentMinActDist* ph->my_type) {
                    coarseDistanceMap->addIntoDistFinal(u,v);
                    toOptimize.push_back(ph);
                }
            } else {
                delete ph;
                host->immaturePoints[i]=0;
            }
        }
    }

    std::vector<PointHessian*> optimized;
    optimized.resize(toOptimize.size());

    if(multiThreading) {
        treadReduce.reduce(
            boost::bind(&FullSystem::activatePointsMT_Reductor, this,
                        &optimized, &toOptimize, _1, _2, _3, _4),
            0, toOptimize.size(), 50
        );
    } else {
        activatePointsMT_Reductor(&optimized, &toOptimize, 0, toOptimize.size(), 0, 0);
    }

    for(unsigned k=0; k<toOptimize.size(); k++) {
        PointHessian* newpoint = optimized[k];
        ImmaturePoint* ph = toOptimize[k];

        if(newpoint != 0 && newpoint != (PointHessian*)((long)(-1))) {
            newpoint->host->immaturePoints[ph->idxInImmaturePoints]=0;
            newpoint->host->pointHessians.push_back(newpoint);
            ef->insertPoint(newpoint);
            for(PointFrameResidual* r : newpoint->residuals) {
                ef->insertResidual(r);
            }
            assert(newpoint->efPoint != 0);
            delete ph;
        }
        else if(newpoint == (PointHessian*)((long)(-1))
                || ph->lastTraceStatus==IPS_OOB)
        {
            delete ph;
            ph->host->immaturePoints[ph->idxInImmaturePoints]=0;
        }
        else
        {
            assert(newpoint == 0 || newpoint == (PointHessian*)((long)(-1)));
        }
    }


    for(FrameHessian* host : frameHessians)
    {
        for(int i=0; i<(int)host->immaturePoints.size(); i++)
        {
            if(host->immaturePoints[i]==0)
            {
                host->immaturePoints[i] = host->immaturePoints.back();
                host->immaturePoints.pop_back();
                i--;
            }
        }
    }
}


void FullSystem::activatePointsOldFirst()
{
    assert(false);
}


void FullSystem::flagPointsForRemoval() {
    assert(EFIndicesValid);

    std::vector<FrameHessian*> fhsToKeepPoints;
    std::vector<FrameHessian*> fhsToMargPoints;

    for(int i=((int)frameHessians.size())-1; i>=0 && i >= ((int)frameHessians.size()); i--) {
        if(!frameHessians[i]->flaggedForMarginalization) {
            fhsToKeepPoints.push_back(frameHessians[i]);
        }
    }

    for(int i=0; i< (int)frameHessians.size(); i++) {
        if(frameHessians[i]->flaggedForMarginalization) {
            fhsToMargPoints.push_back(frameHessians[i]);
        }
    }

    int flag_oob=0, flag_in=0, flag_inin=0, flag_nores=0;

    // go through all active frames
    for(FrameHessian* host : frameHessians) {
        for(unsigned int i=0; i<host->pointHessians.size(); i++) {
            PointHessian* ph = host->pointHessians[i];
            if(ph==0) {
                continue;
            }

            if(ph->idepth < 0 || ph->residuals.size() == 0) {
                host->pointHessiansOut.push_back(ph);
                ph->efPoint->stateFlag = EFPointStatus::PS_DROP;
                host->pointHessians[i]=0;
                flag_nores++;
            } else if(ph->isOOB(fhsToKeepPoints, fhsToMargPoints) ||
                      host->flaggedForMarginalization) {
                flag_oob++;
                if(ph->isInlierNew()) {
                    flag_in++;
                    int ngoodRes = 0;
                    for(PointFrameResidual* r : ph->residuals) {
                        r->resetOOB();
                        r->linearize(&camera_parameters);
                        r->efResidual->isLinearized = false;
                        r->applyRes();
                        if(r->efResidual->isActive()) {
                            r->efResidual->fixLinearizationF(ef);
                            ngoodRes++;
                        }
                    }
                    if(ph->idepth_hessian > setting_minIdepthH_marg) {
                        flag_inin++;
                        ph->efPoint->stateFlag = EFPointStatus::PS_MARGINALIZE;
                        host->pointHessiansMarginalized.push_back(ph);
                    } else {
                        ph->efPoint->stateFlag = EFPointStatus::PS_DROP;
                        host->pointHessiansOut.push_back(ph);
                    }
                } else  {
                    host->pointHessiansOut.push_back(ph);
                    ph->efPoint->stateFlag = EFPointStatus::PS_DROP;
                }

                host->pointHessians[i]=0;
            }
        }

        for(int i=0; i<(int)host->pointHessians.size(); i++) {
            if(host->pointHessians[i] == 0) {
                host->pointHessians[i] = host->pointHessians.back();
                host->pointHessians.pop_back();
                i--;
            }
        }
    }
}


void FullSystem::addActiveFrame(float* image, int id, const float exposure_time) {
    if(isLost) {
        return;
    }

    boost::unique_lock<boost::mutex> lock(trackMutex);

    // =========================== add into allFrameHistory =========================
    // no lock required, as fh is not used anywhere yet.
    FrameShell* shell = new FrameShell(allFrameHistory.size(), id);

    allFrameHistory.push_back(shell);

    // =========================== make Images / derivatives etc. =========================
    FrameHessian* fh = new FrameHessian(image, shell, gamma, exposure_time);

    if(!initialized) {
        // use initializer!
        if(coarseInitializer == 0) {
            // first frame set. fh is kept by coarseInitializer.
            coarseInitializer = new CoarseInitializer(fh, camera_parameters, wG[0], hG[0]);
        } else if(coarseInitializer->trackFrame(fh, outputWrapper))  {
            // if SNAPPED
            // FIXME maybe better to create a new CoarseInitializer instance
            // than editing member variables of the existing one
            initializeFromInitializer(fh);
            lock.unlock();
            deliverTrackedFrame(fh, true);
        }
        return;
    }

    // do front-end operation.
    // =========================== SWAP tracking reference?. =========================
    if(coarseTracker_forNewKF->refFrameID > coarseTracker->refFrameID) {
        boost::unique_lock<boost::mutex> crlock(coarseTrackerSwapMutex);
        CoarseTracker* tmp = coarseTracker;
        coarseTracker = coarseTracker_forNewKF;
        coarseTracker_forNewKF = tmp;
    }

    Vec4 tres = trackNewCoarse(fh);
    if(!std::isfinite((double)tres[0]) || !std::isfinite((double)tres[1]) ||
       !std::isfinite((double)tres[2]) || !std::isfinite((double)tres[3])) {
        printf("Initial Tracking failed: LOST!\n");
        isLost = true;
        return;
    }

    Vec2 refToFh = AffLight::fromToVecExposure(
        coarseTracker->lastRef->ab_exposure,
        fh->ab_exposure,
        coarseTracker->lastRef_aff_g2l, fh->shell->aff_g2l
    );

    // BRIGHTNESS CHECK
    bool needToMakeKF =
       allFrameHistory.size()== 1 ||
       setting_kfGlobalWeight*setting_maxShiftWeightT *  sqrtf((double)tres[1]) /
       (wG[0]+hG[0]) +
       setting_kfGlobalWeight*setting_maxShiftWeightR *  sqrtf((double)tres[2]) /
       (wG[0]+hG[0]) +
       setting_kfGlobalWeight*setting_maxShiftWeightRT * sqrtf((double)tres[3]) /
       (wG[0]+hG[0]) +
       setting_kfGlobalWeight*setting_maxAffineWeight * fabs(logf((float)refToFh[0])) > 1 ||
       2 * coarseTracker->firstCoarseRMSE < tres[0];

    for(IOWrap::Output3DWrapper* ow : outputWrapper) {
        ow->publishCamPose(fh->shell, &camera_parameters);
    }

    lock.unlock();
    deliverTrackedFrame(fh, needToMakeKF);
    return;
}

void FullSystem::deliverTrackedFrame(FrameHessian* fh, bool needKF) {
    if(linearizeOperation) {
        if(goStepByStep && lastRefStopID != coarseTracker->refFrameID)
        {
            MinimalImageF3 img(wG[0], hG[0], fh->dI);
            IOWrap::displayImage("frameToTrack", &img);
            while(true)
            {
                char k = IOWrap::waitKey(0);
                if(k == ' ') break;
            }
            lastRefStopID = coarseTracker->refFrameID;
        }

        if(needKF) {
            makeKeyFrame(fh);
        } else {
            makeNonKeyFrame(fh);
            delete fh;
        }
        return;
    }

    boost::unique_lock<boost::mutex> lock(trackMapSyncMutex);
    unmappedTrackedFrames.push_back(fh);
    if(needKF) {
        needNewKFAfter = fh->shell->trackingRef->id;
    }
    trackedFrameSignal.notify_all();

    while(coarseTracker_forNewKF->refFrameID == -1
          && coarseTracker->refFrameID == -1) {
        mappedFrameSignal.wait(lock);
    }

    lock.unlock();
}

void FullSystem::mappingLoop()
{
    boost::unique_lock<boost::mutex> lock(trackMapSyncMutex);

    while(runMapping) {
        while(unmappedTrackedFrames.size()==0) {
            trackedFrameSignal.wait(lock);
            if(!runMapping) return;
        }

        FrameHessian* fh = unmappedTrackedFrames.front();
        unmappedTrackedFrames.pop_front();


        // guaranteed to make a KF for the very first two tracked frames.
        if(allKeyFramesHistory.size() <= 2) {
            lock.unlock();
            makeKeyFrame(fh);
            lock.lock();
            mappedFrameSignal.notify_all();
            continue;
        }

        if(unmappedTrackedFrames.size() > 3)
            needToKetchupMapping=true;

        // if there are other frames to tracke, do that first.
        if(unmappedTrackedFrames.size() > 0) {
            lock.unlock();
            makeNonKeyFrame(fh);
            delete fh;
            lock.lock();

            if(needToKetchupMapping && unmappedTrackedFrames.size() > 0)
            {
                FrameHessian* fh = unmappedTrackedFrames.front();
                unmappedTrackedFrames.pop_front();
                {
                    // FIXME this part is same as makeNonKeyFrame
                    boost::unique_lock<boost::mutex> crlock(shellPoseMutex);
                    assert(fh->shell->trackingRef != 0);
                    fh->shell->camToWorld = fh->shell->trackingRef->camToWorld *
                                            fh->shell->camToTrackingRef;
                    fh->setEvalPT_scaled(fh->shell->camToWorld.inverse(),fh->shell->aff_g2l);
                }
                delete fh;
            }

        }
        else
        {
            if(setting_realTimeMaxKF || needNewKFAfter >= frameHessians.back()->shell->id)
            {
                lock.unlock();
                makeKeyFrame(fh);
                needToKetchupMapping=false;
                lock.lock();
            }
            else
            {
                lock.unlock();
                makeNonKeyFrame(fh);
                delete fh;
                lock.lock();
            }
        }
        mappedFrameSignal.notify_all();
    }
    printf("MAPPING FINISHED!\n");
}

void FullSystem::blockUntilMappingIsFinished()
{
    boost::unique_lock<boost::mutex> lock(trackMapSyncMutex);
    runMapping = false;
    trackedFrameSignal.notify_all();
    lock.unlock();

    mappingThread.join();
}

void FullSystem::makeNonKeyFrame( FrameHessian* fh) {
    // needs to be set by mapping thread. no lock required since we are in mapping thread.
    boost::unique_lock<boost::mutex> crlock(shellPoseMutex);
    assert(fh->shell->trackingRef != 0);
    fh->shell->camToWorld = fh->shell->trackingRef->camToWorld *
                            fh->shell->camToTrackingRef;
    fh->setEvalPT_scaled(fh->shell->camToWorld.inverse(), fh->shell->aff_g2l);

    traceNewCoarse(fh);
}

bool checkIfInitializationFailed(int nKeyFrames, float rmse) {

    // =========================== Figure Out if INITIALIZATION FAILED =========================
    if(nKeyFrames > 4) {
        return false;
    }
    if(nKeyFrames==2 && rmse > 20*benchmark_initializerSlackFactor) {
        return true;
    }
    if(nKeyFrames==3 && rmse > 13*benchmark_initializerSlackFactor) {
        return true;
    }
    if(nKeyFrames==4 && rmse > 9*benchmark_initializerSlackFactor) {
        return true;
    }
    return false;
}


void FullSystem::makeKeyFrame(FrameHessian* fh) {
    // needs to be set by mapping thread
    makeNonKeyFrame(fh);

    boost::unique_lock<boost::mutex> lock(mapMutex);

    // =========================== Flag Frames to be Marginalized. =========================
    flagFramesForMarginalization();

    // =========================== add New Frame to Hessian Struct. =========================
    fh->idx = frameHessians.size();
    frameHessians.push_back(fh);
    fh->frameID = allKeyFramesHistory.size();
    allKeyFramesHistory.push_back(fh->shell);
    ef->insertFrame(fh);
    ef->setDeltaF(current_camera_parameters);

    setPrecalcValues(frameHessians, createCameraMatrixFromCalibHessian(camera_parameters));

    // =========================== add new residuals for old points =========================
    for(FrameHessian* fh1 : frameHessians) {
        // go through all active frames
        if(fh1 == fh) continue;
        for(PointHessian* ph : fh1->pointHessians) {
            PointFrameResidual* r = new PointFrameResidual(ph, fh1, fh, ResState::IN);
            ph->residuals.push_back(r);
            ef->insertResidual(r);
            ph->lastResiduals[1] = ph->lastResiduals[0];
            ph->lastResiduals[0] = std::pair<PointFrameResidual*, ResState>(r, ResState::IN);
        }
    }

    // =========================== Activate Points (& flag for marginalization). =========================
    activatePointsMT();
    ef->makeIDX();

    // =========================== OPTIMIZE ALL =========================

    fh->frameEnergyTH = frameHessians.back()->frameEnergyTH;
    float rmse = optimize(setting_maxOptIterations);

    // =========================== Figure Out if INITIALIZATION FAILED =========================
    initFailed = checkIfInitializationFailed(allKeyFramesHistory.size(), rmse);
    if(initFailed) {
        printf("I THINK INITIALIZATION FAILED! Resetting.\n");
    }

    if(isLost) return;

    // =========================== REMOVE OUTLIER =========================
    removeOutliers();

    {
        boost::unique_lock<boost::mutex> crlock(coarseTrackerSwapMutex);
        coarseTracker_forNewKF->makeK(&camera_parameters);
        coarseTracker_forNewKF->setCoarseTrackingRef(frameHessians);
    }

    // =========================== (Activate-)Marginalize Points =========================
    flagPointsForRemoval();
    ef->dropPointsF();
    getNullspaces(
        ef->lastNullspaces_pose,
        ef->lastNullspaces_scale);
    ef->marginalizePointsF();

    // =========================== add new Immature points & new residuals =========================
    makeNewTraces(fh, 0);

    for(IOWrap::Output3DWrapper* ow : outputWrapper) {
        ow->publishGraph(ef->connectivityMap);
        ow->publishKeyframes(frameHessians, false, &camera_parameters);
    }

    // =========================== Marginalize Frames =========================

    for(unsigned int i=0; i<frameHessians.size(); i++) {
        if(frameHessians[i]->flaggedForMarginalization) {
            marginalizeFrame(frameHessians[i]);
            i=0;
        }
    }
}


void FullSystem::initializeFromInitializer(FrameHessian* newFrame) {
    boost::unique_lock<boost::mutex> lock(mapMutex);

    // add firstframe.
    FrameHessian* firstFrame = coarseInitializer->firstFrame;

    firstFrame->idx = frameHessians.size();
    frameHessians.push_back(firstFrame);
    firstFrame->frameID = allKeyFramesHistory.size();
    allKeyFramesHistory.push_back(firstFrame->shell);

    ef->insertFrame(firstFrame);
    ef->setDeltaF(current_camera_parameters);

    setPrecalcValues(frameHessians, createCameraMatrixFromCalibHessian(camera_parameters));

    firstFrame->pointHessians.reserve(wG[0] * hG[0] * 0.2f);
    firstFrame->pointHessiansMarginalized.reserve(wG[0] * hG[0] * 0.2f);
    firstFrame->pointHessiansOut.reserve(wG[0] * hG[0] * 0.2f);

    float sumID=1e-5, numID=1e-5;
    for(int i=0; i<coarseInitializer->numPoints[0]; i++) {
        sumID += coarseInitializer->points[0][i].iR;
        numID++;
    }
    float rescaleFactor = 1 / (sumID / numID);

    // randomly sub-select the points I need.
    float keepPercentage = setting_desiredPointDensity / coarseInitializer->numPoints[0];

    if(!setting_debugout_runquiet)
        printf("Initialization: keep %.1f%% (need %d, have %d)!\n", 100*keepPercentage,
               (int)(setting_desiredPointDensity), coarseInitializer->numPoints[0] );

    for(int i=0; i<coarseInitializer->numPoints[0]; i++)
    {
        if(rand()/(float)RAND_MAX > keepPercentage) continue;

        const Pnt point = coarseInitializer->points[0][i];

        ImmaturePoint pt(point.u+0.5f, point.v+0.5f, firstFrame, point.my_type);

        PointHessian* ph = new PointHessian(pt.u, pt.v, pt.my_type, pt.host,
                                            pt.color, pt.weights,
                                            point.iR*rescaleFactor, ph->idepth,
                                            PointHessian::ACTIVE, true);

        firstFrame->pointHessians.push_back(ph);
        ef->insertPoint(ph);
    }

    SE3 firstToNew = coarseInitializer->thisToNext;
    firstToNew.translation() /= rescaleFactor;

    // really no lock required, as we are initializing.
    {
        boost::unique_lock<boost::mutex> crlock(shellPoseMutex);
        firstFrame->setEvalPT_scaled(firstFrame->shell->camToWorld.inverse(),
                                     firstFrame->shell->aff_g2l);
        firstFrame->shell->trackingRef=0;
        firstFrame->shell->camToTrackingRef = SE3();

        newFrame->shell->camToWorld = firstToNew.inverse();

        newFrame->setEvalPT_scaled(newFrame->shell->camToWorld.inverse(),
                                   newFrame->shell->aff_g2l);
        newFrame->shell->trackingRef = firstFrame->shell;
        newFrame->shell->camToTrackingRef = firstToNew.inverse();
    }

    initialized=true;
    printf("INITIALIZE FROM INITIALIZER (%d pts)!\n",
           (int)firstFrame->pointHessians.size());
}

void FullSystem::makeNewTraces(FrameHessian* newFrame, float* gtDepth) {
    int numPointsTotal = pixelSelector->makeMaps(
        newFrame,
        selectionMap,
        setting_desiredImmatureDensity
    );

    newFrame->pointHessians.reserve(numPointsTotal*1.2f);
    //fh->pointHessiansInactive.reserve(numPointsTotal*1.2f);
    newFrame->pointHessiansMarginalized.reserve(numPointsTotal*1.2f);
    newFrame->pointHessiansOut.reserve(numPointsTotal*1.2f);

    for(int y=patternPadding+1; y<hG[0]-patternPadding-2; y++) {
        for(int x=patternPadding+1; x<wG[0]-patternPadding-2; x++) {
            int i = y*wG[0] + x;
            if(selectionMap[i]==0) continue;

            ImmaturePoint* impt;
            try {
                impt = new ImmaturePoint(x, y, newFrame, selectionMap[i]);
            } catch(std::exception &e) {
                std::cerr << e.what() << std::endl;
                continue;
            }

            newFrame->immaturePoints.push_back(impt);
        }
    }
}


void setPrecalcValues(std::vector<FrameHessian*> frameHessians, const Mat33f &K) {
    for(FrameHessian* fh : frameHessians) {
        fh->targetPrecalc.resize(frameHessians.size());
        for(unsigned int i=0; i<frameHessians.size(); i++) {
            fh->targetPrecalc[i].set(fh, frameHessians[i], K);
        }
    }
}

}
