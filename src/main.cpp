// src/main.cpp
#include "render/render.hpp"
#include "scenes/scenes.hpp"
#include <cstdio>
#include <string_view>

// picks the scene from the first argument, ./bsim [plane|spacecraft|sindy|scvx] (plane if there
// isn't one)
int main(int argc, char **argv) {
    const std::string_view scene = argc > 1 ? argv[1] : "plane";
    if (scene != "plane" && scene != "spacecraft" && scene != "sindy" &&
        scene != "scvx") {
        std::fprintf(stderr, "usage: bsim [plane|spacecraft|sindy|scvx]\n");
        return 1;
    }

    BSIM::rendering::init_rendering();
    if (scene == "plane") {
        BSIM::scenes::run_plane_scene();
    } else if (scene == "spacecraft") {
        BSIM::scenes::run_spacecraft_scene();
    } else if (scene == "sindy") {
        BSIM::scenes::run_sindy_scene();
    } else {
        BSIM::scenes::run_scvx_scene();
    }
    CloseWindow();
}
