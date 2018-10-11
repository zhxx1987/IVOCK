set(SOURCE_FILES main.cpp
        read_particles.cpp
        render.cpp
        shadow.cpp
        write_sgi.cpp
        obj_reader.cpp
        bfstream.cpp)

add_IVOCK_executable(renderer ${SOURCE_FILES})