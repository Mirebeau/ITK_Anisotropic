# This is the root ITK CMakeLists file.
cmake_minimum_required(VERSION 2.4)
if(COMMAND CMAKE_POLICY)
  cmake_policy(SET CMP0003 NEW)
endif(COMMAND CMAKE_POLICY)

#set(CMAKE_MACOSX_RPATH 1)


include(UserDirectories.txt)

# This project is designed to be built outside the Insight source tree.
project(MatlabITK)

# Find ITK.
find_package(ITK REQUIRED)
include(${ITK_USE_FILE})

add_definitions(-std=c++11)

add_subdirectory("Headers")
add_subdirectory("../FastMarching/Headers" "${CMAKE_CURRENT_BINARY_DIR}/Useless/FMH")
add_subdirectory("../FastMarching/itkFastMarching_Modified" "${CMAKE_CURRENT_BINARY_DIR}/Useless/FMM")
add_subdirectory("../NormTypes" "${CMAKE_CURRENT_BINARY_DIR}/Useless/NT")
add_subdirectory("../General" "${CMAKE_CURRENT_BINARY_DIR}/Useless/G")
add_subdirectory("../Diffusion/Headers" "${CMAKE_CURRENT_BINARY_DIR}/Useless/D")

include_directories("${MATLAB_APP}${MATLAB_INCLUDE_SUBDIR}")
include_directories(${PROJECT_DIR})
include_directories("Interface")

add_executable(Test Test.cxx ${Matlab_Headers} ${General_Headers} ${NormTypes_Headers} ${AFMIF_Headers})
add_library(MatlabAnisotropicDiffusion SHARED MatlabAnisotropicDiffusion.cxx ${DiffusionOperator_Headers} ${AnisotropicDiff_Headers} ${CoherenceEnhancingDiffusionTensor_Headers})
add_library(MatlabAnisotropicFastMarching_EarlyAbort SHARED MatlabAnisotropicFastMarching_EarlyAbort.cxx ${AnisotropicFM_EarlyAbort_Headers})



source_group("Matlab_Headers" FILES ${Matlab_Headers})
source_group("General_Headers" FILES ${General_Headers})
source_group("NormTypes_Headers" FILES ${NormTypes_Headers})
source_group("AFMIF_Headers" FILES ${AFMIF_Headers})

target_link_libraries(Test ${ITK_LIBRARIES})
target_link_libraries(MatlabAnisotropicDiffusion ${ITK_LIBRARIES})
target_link_libraries(MatlabAnisotropicFastMarching_EarlyAbort ${ITK_LIBRARIES})

set(MATLAB_LIB_DIR "${MATLAB_APP}${MATLAB_LIB_SUBDIR}")
set(MATLAB_LIBRARIES "${MATLAB_LIB_DIR}/libmat.dylib" "${MATLAB_LIB_DIR}/libmex.dylib" "${MATLAB_LIB_DIR}/libmx.dylib")
target_link_libraries(MatlabAnisotropicDiffusion ${MATLAB_LIBRARIES})
target_link_libraries(MatlabAnisotropicFastMarching_EarlyAbort ${MATLAB_LIBRARIES})

file(COPY MFiles DESTINATION "${CMAKE_CURRENT_BINARY_DIR}")
file(COPY "Interface" DESTINATION "${CMAKE_CURRENT_BINARY_DIR}")