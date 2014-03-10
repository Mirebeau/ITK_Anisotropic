# This is the root ITK CMakeLists file.
cmake_minimum_required(VERSION 2.4)
if(COMMAND CMAKE_POLICY)
  cmake_policy(SET CMP0003 NEW)
endif(COMMAND CMAKE_POLICY)


# This project is designed to be built outside the Insight source tree.
project(itkDiffusion)

# Find ITK.
find_package(ITK REQUIRED)
include(${ITK_USE_FILE})

add_subdirectory("Headers")
add_subdirectory("Testing")
add_subdirectory("../NormTypes" "${CMAKE_CURRENT_BINARY_DIR}/NormTypes")
add_subdirectory("../General" "${CMAKE_CURRENT_BINARY_DIR}/General")

set(Project_Headers
   ${NormTypes_Headers}
   ${NormTypes_Implementation}
   ${NormTypes_Base_Headers}
   ${NormTypes_Base_Implementation}
   ${General_Headers}
   ${DiffusionTest_Headers}
   ${DiffusionOperator_Headers}
   ${CoherenceEnhancingDiffusionTensor_Headers}
)

include_directories(${PROJECT_DIR} ${CMAKE_CURRENT_SOURCE_DIR})

add_executable(itkDiffusion itkDiffusion.cxx ${Project_Headers})
add_executable(UserMain "Testing/UserMain.cxx") 

source_group("NormTypes"                 FILES ${NormTypes_Headers})
source_group("NormTypes\\Implementation" FILES ${NormTypes_Implementation})
source_group("NormTypes\\Base"           FILES ${NormTypes_Base_Headers})
source_group("NormTypes\\Base\\Implementation" FILES ${NormTypes_Base_Implementation})

source_group("General" FILES ${General_Headers})


target_link_libraries(itkDiffusion ${ITK_LIBRARIES})
target_link_libraries(UserMain ${ITK_LIBRARIES})

file(COPY "Testing/Images" DESTINATION "${CMAKE_CURRENT_BINARY_DIR}")

