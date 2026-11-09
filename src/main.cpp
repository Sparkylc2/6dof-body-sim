// src/main.cpp
#include "render/render.hpp"
#include "scenes/scenes.hpp"
#include <cstdio>
#include <string_view>

// picks the scene from the first argument, ./bsim [plane|spacecraft] (plane if there isn't one)
int main(int argc, char **argv) {
    const std::string_view scene = argc > 1 ? argv[1] : "plane";
    if (scene != "plane" && scene != "spacecraft") {
        std::fprintf(stderr, "usage: bsim [plane|spacecraft]\n");
        return 1;
    }

    BSIM::rendering::init_rendering();
    if (scene == "plane") {
        BSIM::scenes::run_plane_scene();
    } else {
        BSIM::scenes::run_spacecraft_scene();
    }
    CloseWindow();
}
