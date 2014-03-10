# This is the root ITK CMakeLists file.
cmake_minimum_required(VERSION 2.4)
if(COMMAND CMAKE_POLICY)
  cmake_policy(SET CMP0003 NEW)
endif(COMMAND CMAKE_POLICY)

# This project is designed to be built outside the Insight source tree.
project(ITKFM)

add_definitions(-std=c++11)

#The following lines are required to use Dart
#ENABLE_TESTING()
#INCLUDE(Dart)


# Find ITK.
include(../ITK_DIR.txt)
find_package(ITK REQUIRED)
include(${ITK_USE_FILE})

add_subdirectory("Headers")
add_subdirectory("itkFastMarching_Modified")
add_subdirectory("Testing")
add_subdirectory("../NormTypes" "${CMAKE_CURRENT_BINARY_DIR}/NormTypes")
add_subdirectory("../General" "${CMAKE_CURRENT_BINARY_DIR}/General")

include_directories(${PROJECT_DIR})

set(ITKFM_Headers
   ${NormTypes_Headers}
   ${NormTypes_Implementation}
   ${NormTypes_Base_Headers}
   ${NormTypes_Base_Implementation}
   ${AFMIF_Headers}
   ${AFMIF_Implementation}
   ${itkFastMarching_Modified_Headers}
   ${AFMIF_Tests_Headers}
   ${General_Headers}
)

add_executable(ITKFM ITKFM.cxx ${ITKFM_Headers})
add_executable(UserMain          "Testing/UserMain.cxx")
add_executable(IsotropicTests   "Testing/AFMIFTest.cxx")
add_executable(AnisotropicTests "Testing/AFMIFTest2.cxx" ${AFMIF_Test2_Headers})




source_group("NormTypes"                 FILES ${NormTypes_Headers})
source_group("NormTypes\\Implementation" FILES ${NormTypes_Implementation})
source_group("NormTypes\\Base"           FILES ${NormTypes_Base_Headers})
source_group("NormTypes\\Base\\Implementation" FILES ${NormTypes_Base_Implementation})

source_group("AFMIF\\Implementation"     FILES ${AFMIF_Implementation})
source_group("itkFastMarching_Modified"  FILES ${itkFastMarching_Modified_Headers})
source_group("Tests"                     FILES ${AFMIF_Tests_Headers})

source_group("General" FILES ${General_Headers})
source_group("AFMIF"                     FILES ${AFMIF_Headers})

target_link_libraries(ITKFM ${ITK_LIBRARIES})
target_link_libraries(UserMain ${ITK_LIBRARIES})
target_link_libraries(IsotropicTests ${ITK_LIBRARIES})
target_link_libraries(AnisotropicTests ${ITK_LIBRARIES})


################ Testing ###############
#add_test(HelloWorldTest ITKFM)

