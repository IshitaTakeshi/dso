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

#include "util/NumType.h"
#include "OptimizationBackend/MatrixAccumulators.h"
#include "IOWrapper/Output3DWrapper.h"
#include "util/settings.h"
#include "vector"
#include <math.h>




namespace dso
{
struct CalibHessian;
struct FrameHessian;


struct Pnt
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
    // index in jacobian. never changes (actually, there is no reason why).
    float u,v;

    // idepth / isgood / energy during optimization.
    float idepth;
    bool isGood;
    Vec2f energy;		// (UenergyPhotometric, energyRegularizer)
    bool isGood_new;
    float idepth_new;
    Vec2f energy_new;

    float iR;
    float iRSumNum;

    float lastHessian;
    float lastHessian_new;

    // max stepsize for idepth (corresponding to max. movement in pixel-space).
    float maxstep;

    // idx (x+y*w) of closest point one pyramid level above.
    int parent;
    float parentDist;

    // idx (x+y*w) of up to 10 nearest points in pixel space.
    int neighbours[10];
    float neighboursDist[10];

    float my_type;
    float outlierTH;
};

class CoarseInitializer {
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

    CoarseInitializer(FrameHessian* newFrameHessian,
                      const CalibHessian* HCalib,
                      const int ww_, const int hh_);
    ~CoarseInitializer();

    void setFirst(FrameHessian* newFrameHessian);
    bool trackFrame(FrameHessian* newFrameHessian,
                    std::vector<IOWrap::Output3DWrapper*> &wraps);
    void calcTGrads(FrameHessian* newFrameHessian);

    int frameID;
    bool fixAffine;
    bool printDebug;

    Pnt* points[PYR_LEVELS];
    int numPoints[PYR_LEVELS];
    AffLight thisToNext_aff;
    SE3 thisToNext;

    // TODO make this const
    FrameHessian* firstFrame;
    FrameHessian* newFrame;
private:
    Mat33f Ki[PYR_LEVELS];
    Mat33f K[PYR_LEVELS];
    int w[PYR_LEVELS];
    int h[PYR_LEVELS];
    const int ww; // ugly name
    const int hh; // ugly name
    void makeK(Mat33f K[], Mat33f Ki[],
               const CalibHessian* HCalib, const int w_, const int h_);
    bool snapped;
    int snappedAt;

    // pyramid images & levels on all levels
    Eigen::Vector3f* dINew[PYR_LEVELS];
    Eigen::Vector3f* dIFist[PYR_LEVELS];

    Eigen::DiagonalMatrix<float, 8> wM;

    Vec3f calcResAndGS(
        Vec10f* JbBuffer_new, int lvl,
        Mat88f &H_out, Vec8f &b_out,
        Mat88f &H_out_sc, Vec8f &b_out_sc,
        const SE3 &refToNew, AffLight refToNew_aff,
        float couplingWeight, float alphaW, float alphaK, bool plot);
    Vec3f calcEC(int level, float couplingWeight); // returns OLD NERGY, NEW ENERGY, NUM TERMS.
    void optReg(int lvl, float regWeight = 0.8);

    void propagateUp(int srcLvl);
    void propagateDown(int srcLvl);
    float rescale();

    void resetPoints(int lvl);

    void doStep(Vec10f* JbBuffer, const int level, const float lambda, const Vec8f inc);
    void applyStep(Vec10f* JbBuffer, Vec10f* JbBuffer_new, const int level);

    void makeGradients(Eigen::Vector3f** data);

    void debugPlot(int lvl, std::vector<IOWrap::Output3DWrapper*> &wraps);
    void makeNN();
};




struct FLANNPointcloud
{
    inline FLANNPointcloud() {
        num=0;
        points=0;
    }
    inline FLANNPointcloud(int n, Pnt* p) :  num(n), points(p) {}
    int num;
    Pnt* points;
    inline size_t kdtree_get_point_count() const {
        return num;
    }
    inline float kdtree_distance(const float *p1, const size_t idx_p2,
                                 size_t /*size*/) const
    {
        const float d0=p1[0]-points[idx_p2].u;
        const float d1=p1[1]-points[idx_p2].v;
        return d0*d0+d1*d1;
    }

    inline float kdtree_get_pt(const size_t idx, int dim) const
    {
        if (dim==0) return points[idx].u;
        else return points[idx].v;
    }
    template <class BBOX>
    bool kdtree_get_bbox(BBOX& /* bb */) const {
        return false;
    }
};

}


