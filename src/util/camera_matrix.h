#include "util/NumType.h"
#include "FullSystem/CalibHessian.h"


namespace dso {
Mat33f initializeCameraMatrix(const VecC &camera_parameters);
Mat33f initializeCameraMatrix(const float fx, const float fy,
                              const float cx, const float cy);
Mat33f createCameraMatrixFromCalibHessian(const CalibHessian &HCalib);
}
