set(Matlab_Headers
   ${CMAKE_CURRENT_SOURCE_DIR}/mxIO.h
   ${CMAKE_CURRENT_SOURCE_DIR}/mxIO.hxx
   ${CMAKE_CURRENT_SOURCE_DIR}/mexMessageWrapper.h
   PARENT_SCOPE
)

set(AnisotropicDiff_Headers
 ${CMAKE_CURRENT_SOURCE_DIR}/AnisotropicDiffusion.h
 PARENT_SCOPE
)

set(AnisotropicFM_EarlyAbort_Headers
 ${CMAKE_CURRENT_SOURCE_DIR}/AnisotropicFastMarching_EarlyAbortClass.h
 ${CMAKE_CURRENT_SOURCE_DIR}/AnisotropicFastMarching_EarlyAbort.h
 PARENT_SCOPE
)

set(PROJECT_DIR ${PROJECT_DIR} ${CMAKE_CURRENT_SOURCE_DIR} PARENT_SCOPE)


	
