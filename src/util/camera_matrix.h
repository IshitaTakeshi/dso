#include "util/NumType.h"
#include "FullSystem/CameraParameters.h"


namespace dso {
Mat33f initializeCameraMatrix(const VecC &camera_parameters);
Mat33f initializeCameraMatrix(const float fx, const float fy,
                              const float cx, const float cy);
Mat33f initializeCameraMatrix(const CameraParameters &camera_parameters);
void createCameraMatrixPyramid(Mat33f K[], Mat33f Ki[],
                               const CameraParameters &camera_parameters);
}
