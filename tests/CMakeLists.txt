
include_directories(
    "${CMAKE_SOURCE_DIR}/src"
    "${SC_PATH}/include/common"
)

# trigonometric math with interpolation
add_executable(trig_math_lerp_test trig_math_test.cpp)
target_compile_definitions(trig_math_lerp_test PRIVATE LUT_INTERPOLATION=1)
add_test(NAME trig_math_lerp_test COMMAND trig_math_lerp_test)

# trigonometric math without interpolation
add_executable(trig_math_nolerp_test trig_math_test.cpp)
target_compile_definitions(trig_math_nolerp_test PRIVATE LUT_INTERPOLATION=0)
add_test(NAME trig_math_nolerp_test COMMAND trig_math_nolerp_test)

# window functions
add_executable(window_test window_test.cpp)
add_test(NAME window_test COMMAND window_test)
