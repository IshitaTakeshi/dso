#ifndef UTIL_SCALE_H
#define UTIL_SCALE_H

#include "util/NumType.h"

namespace dso {  //FIXME should not be nested
    const float SCALE_F = 50.0f;
    const float SCALE_C = 50.0f;

    const float SCALE_F_INVERSE = 1.0f / SCALE_F;
    const float SCALE_C_INVERSE = 1.0f / SCALE_C;

    VecC inv_scale_camera_parameters(const VecC &scaled);
    VecC scale_camera_parameters(const VecC &parameters);
}

#endif /* UTIL_SCALE_H */
