set(NormTypes_Base_Headers
   ${CMAKE_CURRENT_SOURCE_DIR}/AdaptiveStencilRefinement2DNormBase.h
   ${CMAKE_CURRENT_SOURCE_DIR}/RiemannianNormBase.h
   PARENT_SCOPE
)

set(NormTypes_Base_Implementation
   ${CMAKE_CURRENT_SOURCE_DIR}/AdaptiveStencilRefinement2DNormBase.hxx
   ${CMAKE_CURRENT_SOURCE_DIR}/RiemannianNormBase.hxx
   PARENT_SCOPE
)

set(PROJECT_DIR ${PROJECT_DIR} ${CMAKE_CURRENT_SOURCE_DIR} PARENT_SCOPE)
