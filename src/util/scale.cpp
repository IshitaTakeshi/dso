#include "util/scale.h"
#include "util/NumType.h"

namespace dso {


// TODO return reference for performance improvement if possible
VecC inv_scale_camera_parameters(const VecC &scaled) {
    VecC parameters;
    parameters[0] = SCALE_F_INVERSE * scaled[0];
    parameters[1] = SCALE_F_INVERSE * scaled[1];
    parameters[2] = SCALE_C_INVERSE * scaled[2];
    parameters[3] = SCALE_C_INVERSE * scaled[3];
    return parameters;
};


// TODO return reference for performance improvement if possible
VecC scale_camera_parameters(const VecC &parameters) {
    VecC scaled;
    scaled[0] = SCALE_F * parameters[0];
    scaled[1] = SCALE_F * parameters[1];
    scaled[2] = SCALE_C * parameters[2];
    scaled[3] = SCALE_C * parameters[3];
    return scaled;
};

}
