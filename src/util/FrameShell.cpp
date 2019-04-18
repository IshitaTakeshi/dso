#include "util/FrameShell.h"

namespace dso {

class FrameShell {
    FrameShell() {
        this->camToWorld = SE3();
        this->aff_g2l = AffLight(0,0);
    }
}
}
