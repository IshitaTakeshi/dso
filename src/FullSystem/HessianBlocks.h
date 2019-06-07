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

#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>

#include "util/NumType.h"
#include "util/MinimalImage.h"
#include "util/gamma.h"
#include "util/scale.h"
#include "FullSystem/CameraParameters.h"
#include "FullSystem/Residuals.h"
#include "FullSystem/ImmaturePoint.h"


namespace dso
{


inline Vec2 affFromTo(const Vec2 &from,
                      const Vec2 &to) {
    // contains affine parameters as XtoWorld.
    return Vec2(from[0] / to[0], (from[1] - to[1]) / to[0]);
}


struct FrameHessian;
struct PointHessian;

class ImmaturePoint;
class FrameShell;

class EFFrame;
class EFPoint;

#define SCALE_XI_ROT 1.0f
#define SCALE_XI_TRANS 0.5f
#define SCALE_F 50.0f
#define SCALE_C 50.0f
#define SCALE_W 1.0f
#define SCALE_A 10.0f
#define SCALE_B 1000.0f

#define SCALE_XI_ROT_INVERSE (1.0f / SCALE_XI_ROT)
#define SCALE_XI_TRANS_INVERSE (1.0f / SCALE_XI_TRANS)
#define SCALE_F_INVERSE (1.0f / SCALE_F)
#define SCALE_C_INVERSE (1.0f / SCALE_C)
#define SCALE_W_INVERSE (1.0f / SCALE_W)
#define SCALE_A_INVERSE (1.0f / SCALE_A)
#define SCALE_B_INVERSE (1.0f / SCALE_B)



struct FrameFramePrecalc {
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
    // static values
    FrameHessian* host;	// defines row
    FrameHessian* target;	// defines column

    // precalc values
    Mat33f PRE_RTll;
    Mat33f PRE_KRKiTll;
    Mat33f PRE_RKiTll;
    Mat33f PRE_RTll_0;

    Vec2f PRE_aff_mode;
    float PRE_b0_mode;

    Vec3f PRE_tTll;
    Vec3f PRE_KtTll;
    Vec3f PRE_tTll_0;

    float distanceLL;

    inline ~FrameFramePrecalc() {}
    inline FrameFramePrecalc() {
        host=target=0;
    }
    void set(FrameHessian* host,
             FrameHessian* target,
             const Mat33f &K);
};


struct FrameHessian {
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
    EFFrame* efFrame;

    // constant info & pre-calculated values
    //DepthImageWrap* frame;
    FrameShell* shell;

    // trace, fine tracking. Used for direction select (not for gradient histograms etc.)
    Eigen::Vector3f* dI;
    // coarse tracking / coarse initializer. NAN in [0] only.
    Eigen::Vector3f* dIp[PYR_LEVELS];
    float* absSquaredGrad[PYR_LEVELS];  // only used for pixel select (histograms etc.). no NAN.

    int frameID;						// incremental ID for keyframes only!
    int idx;

    // Photometric Calibration Stuff
    float frameEnergyTH;	// set dynamically depending on tracking residual
    const float ab_exposure;

    bool flaggedForMarginalization;

    std::vector<PointHessian*> pointHessians;				// contains all ACTIVE points.
    // contains all MARGINALIZED points (= fully marginalized, usually because point went OOB.)
    std::vector<PointHessian*> pointHessiansMarginalized;
    std::vector<PointHessian*> pointHessiansOut;		// contains all OUTLIER points (= discarded.).
    // contains all OUTLIER points (= discarded.).
    std::vector<ImmaturePoint*> immaturePoints;

    Mat66 nullspaces_pose;
    Mat42 nullspaces_affine;
    Vec6 nullspaces_scale;

    // variable info.
    SE3 worldToCam_evalPT;
    Vec10 state_zero;
    Vec10 state_scaled;
    Vec10 state;	// [0-5: worldToCam-leftEps. 6-7: a,b]
    Vec10 step;
    Vec10 step_backup;


    EIGEN_STRONG_INLINE const SE3 &get_worldToCam_evalPT() const {
        return worldToCam_evalPT;
    }
    EIGEN_STRONG_INLINE const Vec10 &get_state_zero() const {
        return state_zero;
    }
    EIGEN_STRONG_INLINE const Vec10 &get_state() const {
        return state;
    }
    EIGEN_STRONG_INLINE const Vec10 &get_state_scaled() const {
        return state_scaled;
    }
    EIGEN_STRONG_INLINE const Vec10 get_state_minus_stateZero() const {
        return get_state() - get_state_zero();
    }


    // precalc values
    SE3 PRE_worldToCam;
    SE3 PRE_camToWorld;
    std::vector<
        FrameFramePrecalc,
        Eigen::aligned_allocator<FrameFramePrecalc>
        > targetPrecalc;
    MinimalImageB3* debugImage;

    inline Vec6 w2c_leftEps() const {
        return get_state_scaled().head<6>();
    }
    inline AffLight aff_g2l() const {
        return AffLight(get_state_scaled()[6], get_state_scaled()[7]);
    }
    inline AffLight aff_g2l_0() const {
        return AffLight(get_state_zero()[6]*SCALE_A, get_state_zero()[7]*SCALE_B);
    }

    void setStateZero(const Vec10 &state_zero);
    void setState(const Vec10 &state);

    inline void setEvalPT(const SE3 &worldToCam_evalPT, const Vec10 &state)
    {

        this->worldToCam_evalPT = worldToCam_evalPT;
        setState(state);
        setStateZero(state);
    };

    inline void setEvalPT_scaled(const SE3 &worldToCam_evalPT,
                                 const AffLight &aff_g2l) {
        Vec10 initial_state = Vec10::Zero();
        initial_state[6] = SCALE_A_INVERSE * aff_g2l.a;
        initial_state[7] = SCALE_B_INVERSE * aff_g2l.b;

        this->worldToCam_evalPT = worldToCam_evalPT;
        setState(initial_state);
        setStateZero(this->get_state());
    };

    void release();

    ~FrameHessian();
    FrameHessian(float* image, FrameShell* shell_,
                 const Gamma &gamma, const float ab_exposure_);

    void makeImages(const float* color, const Gamma &gamma);

    inline Vec10 getPrior()
    {
        Vec10 p =  Vec10::Zero();
        if(frameID==0)
        {
            p.head<3>() = Vec3::Constant(setting_initialTransPrior);
            p.segment<3>(3) = Vec3::Constant(setting_initialRotPrior);
            if(setting_solverMode & SOLVER_REMOVE_POSEPRIOR) p.head<6>().setZero();

            p[6] = setting_initialAffAPrior;
            p[7] = setting_initialAffBPrior;
        }
        else
        {
            if(setting_affineOptModeA < 0)
                p[6] = setting_initialAffAPrior;
            else
                p[6] = setting_affineOptModeA;

            if(setting_affineOptModeB < 0)
                p[7] = setting_initialAffBPrior;
            else
                p[7] = setting_affineOptModeB;
        }
        p[8] = setting_initialAffAPrior;
        p[9] = setting_initialAffBPrior;
        return p;
    }

    inline Vec10 getPriorZero()
    {
        return Vec10::Zero();
    }

};


// hessian component associated with one point.
struct PointHessian {
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
    EFPoint* efPoint;

    // static values
    float color[MAX_RES_PER_POINT];			// colors in host frame
    float weights[MAX_RES_PER_POINT];		// host-weights for respective residuals.

    float u,v;
    int idx;
    float energyTH; // TODO make this const
    FrameHessian* host;
    bool hasDepthPrior;

    float my_type;

    float idepth_zero;
    float idepth;
    float step;
    float step_backup;

    float idepth_hessian;
    float maxRelBaseline;
    int numGoodResiduals;

    enum PtStatus {ACTIVE=0, INACTIVE, OUTLIER, OOB, MARGINALIZED};
    PtStatus status;

    inline void setPointStatus(PtStatus s) {
        status=s;
    }

    inline void setIdepth(float idepth) {
        this->idepth = idepth;
    }

    inline void setIdepthZero(float idepth) {
        idepth_zero = idepth;
    }

    // only contains good residuals (not OOB and not OUTLIER). Arbitrary order.
    std::vector<PointFrameResidual*> residuals;
    // contains information about residuals to the last two (!) frames.
    // ([0] = latest, [1] = the one before).
    std::pair<PointFrameResidual*, ResState> lastResiduals[2];

    void release();

    PointHessian(const float u_, const float v_,
                 const float my_type_, FrameHessian* host_,
                 const float color_[], const float weights_[],
                 const float idepth_, const float idepth_zero_,
                 const PtStatus status_, const bool hasDepthPrior_ = false);

    inline ~PointHessian() {
        assert(efPoint==0);
        release();
    }


    inline bool isOOB(const std::vector<FrameHessian*>& toKeep,
                      const std::vector<FrameHessian*>& toMarg) const
    {

        int visInToMarg = 0;
        for(PointFrameResidual* r : residuals) {
            if(r->state_state != ResState::IN) {
                continue;
            }
            for(FrameHessian* k : toMarg) {
                if(r->target == k) {
                    visInToMarg++;
                }
            }
        }
        if((int)residuals.size() >= setting_minGoodActiveResForMarg &&
           numGoodResiduals > setting_minGoodResForMarg+10 &&
           (int)residuals.size()-visInToMarg < setting_minGoodActiveResForMarg) {
            return true;
        }

        if(lastResiduals[0].second == ResState::OOB) {
            return true;
        }

        if(residuals.size() < 2) {
            return false;
        }

        if(lastResiduals[0].second == ResState::OUTLIER &&
           lastResiduals[1].second == ResState::OUTLIER) {
            return true;
        }

        return false;
    }


    inline bool isInlierNew()
    {
        return (int)residuals.size() >= setting_minGoodActiveResForMarg
               && numGoodResiduals >= setting_minGoodResForMarg;
    }

};





}

