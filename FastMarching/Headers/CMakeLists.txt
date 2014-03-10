set(AFMIF_Headers 
   ${CMAKE_CURRENT_SOURCE_DIR}/itkAnisotropicFastMarchingImageFilter.h
   ${CMAKE_CURRENT_SOURCE_DIR}/CommandLineCall2.h
   PARENT_SCOPE
)

set(AFMIF_Implementation 
   ${CMAKE_CURRENT_SOURCE_DIR}/itkAnisotropicFastMarchingImageFilter.hxx
   ${CMAKE_CURRENT_SOURCE_DIR}/itkAnisotropicFastMarchingImageFilter_StencilData.hxx
   ${CMAKE_CURRENT_SOURCE_DIR}/itkAnisotropicFastMarchingImageFilter_Geodesic.hxx
   PARENT_SCOPE
)

set(PROJECT_DIR ${PROJECT_DIR} ${CMAKE_CURRENT_SOURCE_DIR} PARENT_SCOPE)

