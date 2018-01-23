add_subdirectory(Base)

set(NormTypes_Base_Headers ${NormTypes_Base_Headers} PARENT_SCOPE)
set(NormTypes_Base_Implementation ${NormTypes_Base_Implementation} PARENT_SCOPE)


set(NormTypes_Headers
   ${CMAKE_CURRENT_SOURCE_DIR}/Riemannian2DNorm.h
   ${CMAKE_CURRENT_SOURCE_DIR}/Riemannian3DNorm.h
   ${CMAKE_CURRENT_SOURCE_DIR}/Finsler2DNorm.h
   ${CMAKE_CURRENT_SOURCE_DIR}/AsymmetricQuadratic2DNorm.h
   ${CMAKE_CURRENT_SOURCE_DIR}/ExtendedNorm.h
   PARENT_SCOPE
)

set(NormTypes_Implementation
   ${CMAKE_CURRENT_SOURCE_DIR}/Riemannian2DNorm.hxx
   ${CMAKE_CURRENT_SOURCE_DIR}/Riemannian3DNorm.hxx
   ${CMAKE_CURRENT_SOURCE_DIR}/Finsler2DNorm.hxx
   ${CMAKE_CURRENT_SOURCE_DIR}/AsymmetricQuadratic2DNorm.hxx
   ${CMAKE_CURRENT_SOURCE_DIR}/ExtendedNorm.hxx
   PARENT_SCOPE
)


set(PROJECT_DIR ${PROJECT_DIR} ${CMAKE_CURRENT_SOURCE_DIR} PARENT_SCOPE)