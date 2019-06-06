
#ifndef FULLSYSTEM_CALIBHESSIAN_H
#define FULLSYSTEM_CALIBHESSIAN_H

#include "util/NumType.h"

namespace dso {
struct CameraParameters {
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

    VecC camera_parameters;

    inline ~CameraParameters() {
    }

    inline CameraParameters(const Eigen::Matrix3f &K) {
        this->camera_parameters << K(0,0), K(1,1), K(0,2), K(1,2);
    };

    VecC get() const {
        return this->camera_parameters;
    }

    // normal mode: use the optimized parameters everywhere!
    inline float fxl() const {
        return (float)(this->camera_parameters[0]);
    }
    inline float fyl() const {
        return (float)(this->camera_parameters[1]);
    }
    inline float cxl() const {
        return (float)(this->camera_parameters[2]);
    }
    inline float cyl() const {
        return (float)(this->camera_parameters[3]);
    }
    inline float const fxli() const {
        return 1.0f / fxl();
    }
    inline float const fyli() const {
        return 1.0f / fyl();
    }
    inline float const cxli() const {
        return -cxl() / fxl();
    }
    inline float const cyli() const {
        return -cyl() / fyl();
    }
};

}

#endif /* FULLSYSTEM_CALIBHESSIAN_H */
