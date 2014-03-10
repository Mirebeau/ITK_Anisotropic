set(DiffusionOperator_Headers
   ${CMAKE_CURRENT_SOURCE_DIR}/AnisotropicDiffusionLBRImageFilter.h
   ${CMAKE_CURRENT_SOURCE_DIR}/AnisotropicDiffusionLBRImageFilter.hxx
   ${CMAKE_CURRENT_SOURCE_DIR}/LinearAnisotropicDiffusionLBRImageFilter.h
   ${CMAKE_CURRENT_SOURCE_DIR}/LinearAnisotropicDiffusionLBRImageFilter.hxx
   PARENT_SCOPE
)

set(CoherenceEnhancingDiffusionTensor_Headers 
   ${CMAKE_CURRENT_SOURCE_DIR}/CoherenceEnhancingDiffusionFilter.h
   PARENT_SCOPE
)

set(PROJECT_DIR ${PROJECT_DIR} ${CMAKE_CURRENT_SOURCE_DIR} PARENT_SCOPE)

