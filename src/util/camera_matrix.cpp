#include "util/camera_matrix.h"
#include "util/NumType.h"


namespace dso {

Mat33f initializeCameraMatrix(const float fx, const float fy,
                              const float cx, const float cy) {
    Mat33f K = Mat33f::Zero();
    K(0,0) = fx;
    K(1,1) = fy;
    K(0,2) = cx;
    K(1,2) = cy;
    K(2,2) = 1;
    return K;
}


Mat33f createCameraMatrixFromCalibHessian(const CalibHessian &HCalib) {
    return initializeCameraMatrix(HCalib.fxl(), HCalib.fyl(),
                                  HCalib.cxl(), HCalib.cyl());
}

}
