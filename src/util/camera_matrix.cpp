#include "util/camera_matrix.h"
#include "util/NumType.h"
#include "util/settings.h"


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


Mat33f initializeCameraMatrix(const CameraParameters &camera_parameters) {
    return initializeCameraMatrix(camera_parameters.fxl(), camera_parameters.fyl(),
                                  camera_parameters.cxl(), camera_parameters.cyl());
}


void createCameraMatrixPyramid(Mat33f K[], Mat33f Ki[],
                               const CameraParameters &camera_parameters) {

    for (int level = 0; level < pyrLevelsUsed; ++ level) {
        float L = std::pow(2, level);
        float fx = camera_parameters.fxl() / L;
        float fy = camera_parameters.fyl() / L;
        float cx = (camera_parameters.cxl() + 0.5) / L - 0.5;
        float cy = (camera_parameters.cyl() + 0.5) / L - 0.5;

        K[level] = initializeCameraMatrix(fx, fy, cx, cy);
        Ki[level] = K[level].inverse();
    }
}

}
