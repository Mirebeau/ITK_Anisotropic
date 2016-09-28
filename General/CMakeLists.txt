set(General_Headers
   ${CMAKE_CURRENT_SOURCE_DIR}/IsoFM_Geodesic.h
   ${CMAKE_CURRENT_SOURCE_DIR}/UnrollingConstructors.h
   ${CMAKE_CURRENT_SOURCE_DIR}/ExportToMathematica.h
   ${CMAKE_CURRENT_SOURCE_DIR}/RescaleAndExportBMP.h
   ${CMAKE_CURRENT_SOURCE_DIR}/RescaleAndExportBMP.hxx
   ${CMAKE_CURRENT_SOURCE_DIR}/StructureTensorImageFilter.h
   ${CMAKE_CURRENT_SOURCE_DIR}/SyntheticImageSource.h
   ${CMAKE_CURRENT_SOURCE_DIR}/UnaryFunctorWithIndexImageFilter.h
   ${CMAKE_CURRENT_SOURCE_DIR}/Macro.h
   PARENT_SCOPE
)

set(PROJECT_DIR ${PROJECT_DIR} ${CMAKE_CURRENT_SOURCE_DIR} PARENT_SCOPE)
