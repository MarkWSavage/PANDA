cmake_minimum_required(VERSION 3.16 FATAL_ERROR)
project(PANDA)

# Find Geant4 package
find_package(Geant4 REQUIRED ui_all vis_all)

# Include Geant4 settings
include(${Geant4_USE_FILE})

# Include directories
include_directories(${PROJECT_SOURCE_DIR}/include)

# Source files
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

# Build executable
add_executable(PANDA PANDA.cc ${sources} ${headers})

# Link Geant4 libraries
target_link_libraries(PANDA ${Geant4_LIBRARIES})

# Copy macro files to build directory
file(GLOB MACRO_FILES "*.mac")
file(COPY ${MACRO_FILES} DESTINATION ${PROJECT_BINARY_DIR})

# Install target
install(TARGETS PANDA DESTINATION bin)
