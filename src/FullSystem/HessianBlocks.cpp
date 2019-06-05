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



#include "FullSystem/HessianBlocks.h"
#include "util/FrameShell.h"
#include "util/math.h"


namespace dso {


// FIXME passing the reference of 'host'. Dengerous
PointHessian::PointHessian(const float u_, const float v_,
                           const float my_type_, FrameHessian* host_,
                           const float color_[], const float weights_[],
                           const float idepth_min, const float idepth_max) :
    u(u_), v(v_), my_type(my_type_), host(host_) {

    assert(std::isfinite(idepth_max));
    hasDepthPrior=false;

    idepth_hessian=0;
    maxRelBaseline=0;
    numGoodResiduals=0;

    setIdepth((idepth_max + idepth_min)*0.5);
    setPointStatus(PointHessian::INACTIVE);

    memcpy(color, color_, sizeof(float)*patternNum);
    memcpy(weights, weights_, sizeof(float)*patternNum);

    efPoint=0;
}


void PointHessian::release()
{
    for(unsigned int i=0; i<residuals.size(); i++) delete residuals[i];
    residuals.clear();
}

void FrameHessian::setState(const Vec10 &state) {
    this->state = state;
    state_scaled.segment<3>(0) = SCALE_XI_TRANS * state.segment<3>(0);
    state_scaled.segment<3>(3) = SCALE_XI_ROT * state.segment<3>(3);
    state_scaled[6] = SCALE_A * state[6];
    state_scaled[7] = SCALE_B * state[7];
    state_scaled[8] = SCALE_A * state[8];
    state_scaled[9] = SCALE_B * state[9];

    PRE_worldToCam = SE3::exp(w2c_leftEps()) * get_worldToCam_evalPT();
    PRE_camToWorld = PRE_worldToCam.inverse();
    //setCurrentNullspace();
};

void FrameHessian::setStateZero(const Vec10 &state_zero) {
    assert(state_zero.head<6>().squaredNorm() < 1e-20);

    this->state_zero = state_zero;


    for(int i=0; i<6; i++) {
        Vec6 eps;
        eps.setZero();
        eps[i] = 1e-3;
        SE3 EepsP = Sophus::SE3::exp(eps);
        SE3 EepsM = Sophus::SE3::exp(-eps);
        SE3 w2c_leftEps_P_x0 = (get_worldToCam_evalPT() * EepsP) *
                               get_worldToCam_evalPT().inverse();
        SE3 w2c_leftEps_M_x0 = (get_worldToCam_evalPT() * EepsM) *
                               get_worldToCam_evalPT().inverse();
        nullspaces_pose.col(i) = (w2c_leftEps_P_x0.log() - w2c_leftEps_M_x0.log())/
                                 (2e-3);
    }
    //nullspaces_pose.topRows<3>() *= SCALE_XI_TRANS_INVERSE;
    //nullspaces_pose.bottomRows<3>() *= SCALE_XI_ROT_INVERSE;

    // scale change
    SE3 w2c_leftEps_P_x0 = (get_worldToCam_evalPT());
    w2c_leftEps_P_x0.translation() *= 1.00001;
    w2c_leftEps_P_x0 = w2c_leftEps_P_x0 * get_worldToCam_evalPT().inverse();
    SE3 w2c_leftEps_M_x0 = (get_worldToCam_evalPT());
    w2c_leftEps_M_x0.translation() /= 1.00001;
    w2c_leftEps_M_x0 = w2c_leftEps_M_x0 * get_worldToCam_evalPT().inverse();
    nullspaces_scale = (w2c_leftEps_P_x0.log() - w2c_leftEps_M_x0.log())/(2e-3);


    nullspaces_affine.setZero();
    nullspaces_affine.topLeftCorner<2,1>()  = Vec2(1,0);
    assert(ab_exposure > 0);
    nullspaces_affine.topRightCorner<2,1>() = Vec2(0, expf(aff_g2l_0().a)*ab_exposure);
};



void FrameHessian::release() {
    // DELETE POINT
    // DELETE RESIDUAL
    for(unsigned int i=0; i<pointHessians.size(); i++) {
        delete pointHessians[i];
    }
    for(unsigned int i=0; i<pointHessiansMarginalized.size(); i++) {
        delete pointHessiansMarginalized[i];
    }
    for(unsigned int i=0; i<pointHessiansOut.size(); i++) {
        delete pointHessiansOut[i];
    }
    for(unsigned int i=0; i<immaturePoints.size(); i++) {
        delete immaturePoints[i];
    }

    pointHessians.clear();
    pointHessiansMarginalized.clear();
    pointHessiansOut.clear();
    immaturePoints.clear();
}

FrameHessian::~FrameHessian()
{
    assert(efFrame==0);
    release();
    for(int i=0; i<pyrLevelsUsed; i++)
    {
        delete[] dIp[i];
        delete[]  absSquaredGrad[i];

    }

    if(debugImage != 0) delete debugImage;
};


// TODO swap names ab_exposure_  <--> ab_exposure
FrameHessian::FrameHessian(float* image, FrameShell* shell_,
                           const Gamma &gamma, const float ab_exposure_) :
    shell(shell_), ab_exposure(ab_exposure_) {
    flaggedForMarginalization=false;
    frameID = -1;
    efFrame = 0;
    frameEnergyTH = 8*8*patternNum;

    debugImage=0;

    this->makeImages(image, gamma);
};


void FrameHessian::makeImages(const float* color, const Gamma &gamma) {
    for(int i=0; i<pyrLevelsUsed; i++) {
        dIp[i] = new Eigen::Vector3f[wG[i]*hG[i]];
        absSquaredGrad[i] = new float[wG[i]*hG[i]];
    }
    dI = dIp[0];

    // make d0
    int w = wG[0];
    int h = hG[0];
    for(int i=0; i<w*h; i++) {
        dI[i][0] = color[i];
    }

    for(int level=0; level<pyrLevelsUsed; level++) {
        int wl = wG[level], hl = hG[level];
        Eigen::Vector3f* dI_l = dIp[level];

        if(level>0) {
            int wlm1 = wG[level-1];
            Eigen::Vector3f* dI_lm = dIp[level-1];

            for(int y=0; y<hl; y++) {
                for(int x=0; x<wl; x++) {
                    dI_l[x + y*wl][0] = 0.25f * (
                        dI_lm[2*x   + 2*y*wlm1][0] +
                        dI_lm[2*x+1 + 2*y*wlm1][0] +
                        dI_lm[2*x   + 2*y*wlm1+wlm1][0] +
                        dI_lm[2*x+1 + 2*y*wlm1+wlm1][0]
                    );
                }
            }
        }

        float* dabs_l = absSquaredGrad[level];

        for(int idx=wl; idx < wl*(hl-1); idx++) {
            float dx = 0.5f*(dI_l[idx+1][0] - dI_l[idx-1][0]);
            float dy = 0.5f*(dI_l[idx+wl][0] - dI_l[idx-wl][0]);

            if(!std::isfinite(dx)) dx=0;
            if(!std::isfinite(dy)) dy=0;

            dI_l[idx][1] = dx;
            dI_l[idx][2] = dy;

            dabs_l[idx] = dx*dx+dy*dy;

            if(setting_gammaWeightsPixelSelect==1 && gamma.isAvailable) {
                int c = (int)(clamp(std::round((float)(dI_l[idx][0])), 5.0, 250.0));
                float gw = gamma.get(c+1) - gamma.get(c);
                // convert to gradient of original color space (before removing response).
                dabs_l[idx] *= gw*gw;
            }
        }
    }
}

void FrameFramePrecalc::set(FrameHessian* host, FrameHessian* target,
                            const Mat33f &K) {
    this->host = host;
    this->target = target;

    SE3 leftToLeft_0 = target->get_worldToCam_evalPT() *
                       host->get_worldToCam_evalPT().inverse();
    PRE_RTll_0 = (leftToLeft_0.rotationMatrix()).cast<float>();
    PRE_tTll_0 = (leftToLeft_0.translation()).cast<float>();

    SE3 leftToLeft = target->PRE_worldToCam * host->PRE_camToWorld;
    PRE_RTll = (leftToLeft.rotationMatrix()).cast<float>();
    PRE_tTll = (leftToLeft.translation()).cast<float>();
    distanceLL = leftToLeft.translation().norm();

    PRE_KRKiTll = K * PRE_RTll * K.inverse();
    PRE_RKiTll = PRE_RTll * K.inverse();
    PRE_KtTll = K * PRE_tTll;

    PRE_aff_mode = AffLight::fromToVecExposure(
        host->ab_exposure,
        target->ab_exposure,
        host->aff_g2l(),
        target->aff_g2l()).cast<float>();
    PRE_b0_mode = host->aff_g2l_0().b;
}

}
