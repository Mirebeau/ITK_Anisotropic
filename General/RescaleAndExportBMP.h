//
//  RescaleAndExportBMP.h
//  itkDiffusion
//
//  Created by MIREBEAU Jean-Marie on 20/11/12.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#ifndef itkDiffusion_RescaleAndExportBMP_h
#define itkDiffusion_RescaleAndExportBMP_h

//#include "itkImageToImageFilter.h"

#include "itkLightObject.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkImageFileWriter.h" 


// With an itk class ...

namespace itk {
    
    template<typename TImage> 
    class RescaleAndExportBMP : public LightObject
    {
    public:
        typedef RescaleAndExportBMP Self; 
        typedef LightObject Superclass; 
        typedef SmartPointer<Self> Pointer; 
        typedef SmartPointer<const Self> ConstPointer;
        
        // Method for creation through the object factory. 
        itkNewMacro(Self);
        // Run-time type information (and related methods).  
        itkTypeMacro(RescaleAndExportBMP, LightObject);
        
        typedef TImage ImageType; 
        
        void SetFileName(const std::string & filename){m_filename = filename;}
        
        static const unsigned int Dimension = TImage::ImageDimension;
        
        RescaleAndExportBMP();
        
        void Update();
        void SetInput(const ImageType *image){m_caster->SetInput(image);};
    protected:
        std::string m_filename;
        
        typedef unsigned char                            OutputPixelType;
        typedef itk::Image< OutputPixelType, Dimension > OutputImageType;
        
        typedef itk::RescaleIntensityImageFilter<ImageType, OutputImageType >   CastFilterType;
        typename CastFilterType::Pointer m_caster;
        
        typedef  itk::ImageFileWriter< OutputImageType > WriterType;
        typename WriterType::Pointer m_writer;
    };
    
}

#ifndef ITK_MANUAL_INSTANTIATION
#include "RescaleAndExportBMP.hxx"
#endif

#endif
