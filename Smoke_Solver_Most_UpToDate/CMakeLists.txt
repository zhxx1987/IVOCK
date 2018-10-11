find_package(TBB REQUIRED)

set(SOURCE_FILES main.cpp
        Smoke_solver3D.cpp
        GeometricLevelGen.cpp)

add_IVOCK_executable(IVOCK_SIM ${SOURCE_FILES})
target_link_libraries(IVOCK_SIM PRIVATE ${TBB_LIBRARIES})
