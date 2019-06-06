#include "util/camera_matrix.h"
#include "util/NumType.h"


namespace dso {

Mat33f initializeCameraMatrix(const VecC &camera_parameters) {
    return initializeCameraMatrix(
        camera_parameters[0],
        camera_parameters[1],
        camera_parameters[2],
        camera_parameters[3]
    );
}

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


Mat33f createCameraMatrixFromCalibHessian(const CameraParameters &camera_parameters) {
    return initializeCameraMatrix(camera_parameters.fxl(), camera_parameters.fyl(),
                                  camera_parameters.cxl(), camera_parameters.cyl());
}

}
