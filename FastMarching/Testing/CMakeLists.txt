set(AFMIF_Tests_Headers 
   ${CMAKE_CURRENT_SOURCE_DIR}/NonRegressionTests.h  
   PARENT_SCOPE
)

set(AFMIF_Tests2_Headers 
	${CMAKE_CURRENT_SOURCE_DIR}/Tutorial.h  
	${CMAKE_CURRENT_SOURCE_DIR}/Tutorial3D.h  
	${CMAKE_CURRENT_SOURCE_DIR}/TubularSegmentation2DTutorial.h 
	${CMAKE_CURRENT_SOURCE_DIR}/CurveNeighborhoodSource.h.h  
	${CMAKE_CURRENT_SOURCE_DIR}/CurveNeighborhoodSource.h	${CMAKE_CURRENT_SOURCE_DIR}/CurveNeighborhoodSource.hxx	${CMAKE_CURRENT_SOURCE_DIR}/TubeCenterline2D_Image.h	${CMAKE_CURRENT_SOURCE_DIR}/TubeCenterline2D_Test.h	${CMAKE_CURRENT_SOURCE_DIR}/TubeCenterline3D_Image.h	${CMAKE_CURRENT_SOURCE_DIR}/TubeCenterline3D_Test.h	${CMAKE_CURRENT_SOURCE_DIR}/TubularBand_Test.h
	${CMAKE_CURRENT_SOURCE_DIR}/GenerateTestImages.h	${CMAKE_CURRENT_SOURCE_DIR}/Segmentation_Image.h	${CMAKE_CURRENT_SOURCE_DIR}/Segmentation_Test.h	${CMAKE_CURRENT_SOURCE_DIR}/Topographic_Image.h	${CMAKE_CURRENT_SOURCE_DIR}/Topographic_Test.h	${CMAKE_CURRENT_SOURCE_DIR}/Tubular3D_Test.h
   PARENT_SCOPE
)


set(PROJECT_DIR ${PROJECT_DIR} ${CMAKE_CURRENT_SOURCE_DIR} PARENT_SCOPE)
