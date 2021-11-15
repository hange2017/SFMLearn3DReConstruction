project(class2)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "-fPIC")

include_directories("../..")

#test math basic
set(MATH_BASIC_FILE
        task1-1_test_math_basic.cc)
add_executable(task1_test_math_basic ${MATH_BASIC_FILE} )
target_link_libraries(task1_test_math_basic  util core )

#home work task1-2_test_camera
set(task1_2 task1-2_test_camera.cc)
add_executable(task1_2 ${task1_2})
target_link_libraries(task1_2 util core)


# test camera model
set(CAMERA_MODEL
        task1-2_test_camera.cc)

add_executable(task2_test_camera_model ${CAMERA_MODEL}  )
target_link_libraries(task2_test_camera_model util core features )


# test feature matching
set(FEATURE_MATCHING
        task1-6_test_matching.cc
        visualizer.h
        visualizer.cc)

add_executable(task3_test_feature_matching ${FEATURE_MATCHING})
target_link_libraries(task3_test_feature_matching sfm util core features  )


# test fundamental 8 point
set(FUNDAMENTAL_8_POINT_FILE
        task1-3_test_fundamental_8_point.cc)
add_executable(task4_test_fundamental_8_point ${FUNDAMENTAL_8_POINT_FILE} )
target_link_libraries(task4_test_fundamental_8_point sfm util core features )


# test fundamental ransac
set(FUNDAMENTAL_RANSAC_FILE
        task1-4_test_fundamental_ransac.cc)
add_executable(task5_test_fundamental_ransac ${FUNDAMENTAL_RANSAC_FILE})
target_link_libraries(task5_test_fundamental_ransac sfm util core features )


set(POSE_FROM_FUNDAMENTAL
        task1-5_test_pose_from_fundamental.cc
        )
add_executable(task6_test_pose_from_fundamental ${POSE_FROM_FUNDAMENTAL} )
target_link_libraries(task6_test_pose_from_fundamental sfm util core features )

set(TEST_FEATURES task1-7_test_features.cc )
add_executable(task7_test_features ${TEST_FEATURES})
target_link_libraries(task7_test_features sfm util core features)