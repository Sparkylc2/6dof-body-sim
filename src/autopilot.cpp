#include "autopilot.hpp"
#include "force_models.hpp"
#include "frames.hpp"
#include "sensors.hpp" // ground_speed_course

namespace BSIM {

ControlState true_control_state(const AircraftState &x, const AirState &air) {
    ControlState s;

    s.pn = x.pos().x();
    s.pe = x.pos().y();
    s.h = -x.pos().z();

    const AirData air_info = air_data(x.vel(), x.quat(), air);
    s.Va = air_info.Va;
    s.alpha = air_info.alpha;
    s.beta = air_info.beta;

    const Vector3d euler = quat_to_euler(x.quat());
    s.phi = euler.x();
    s.theta = euler.y();
    s.psi = euler.z();

    // course is where the ground velocity points, so it differs by the wind and sideslip
    const Vector2d Vg_chi = ground_speed_course(rotation(x.quat()) * x.vel());
    s.Vg = Vg_chi.x();
    s.chi = Vg_chi.y();

    s.p = x.rates().x();
    s.q = x.rates().y();
    s.r = x.rates().z();

    return s;
}

const FixedWingControls &Autopilot::step(double t, const AutopilotCommands &cmd,
                                         const ControlState &x) {
    if (clock.ready(t)) {
        controls = update(cmd, x);
    }
    return controls;
}

} // namespace BSIM
