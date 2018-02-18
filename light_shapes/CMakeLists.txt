cmake_minimum_required(VERSION 3.9)
project(light_shapes)

set(CMAKE_CXX_STANDARD 11)

add_executable(light_shapes
        3d_clip.c
        backface_elimination.c
        D3d_matrix.c
        D3d_matrix.h
        D3d_matrixS.c
        example_painters.c
        lighting.c
        lightlab.c
        lighty_boi.c
        lit_yellow_sphere.c
        scene01a_nu01.c
        scene01a_nu03.c
        scene_kid.c
        scene_kid_done.c
        scene_works.c
        sphere.c
        sphere_kid.c
        stickman_with_matricesA.c
        stickman_with_matricesB.c
        view_test04.c
        view_test05.c)
