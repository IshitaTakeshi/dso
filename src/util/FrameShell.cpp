#include "util/FrameShell.h"

namespace dso {

FrameShell::FrameShell(int marginalizedAt_, int incoming_id_) :
        marginalizedAt(marginalizedAt_), id(marginalizedAt_), incoming_id(incoming_id_),
        camToWorld(SE3()), camToTrackingRef(SE3()), aff_g2l(AffLight(0,0)), poseValid(true),
        movedByOpt(0), trackingRef(0) {
}
}
