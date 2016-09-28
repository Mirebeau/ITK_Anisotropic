//
//  RescaleAndExport_Function.h
//  ITKFM
//
//  Created by Jean-Marie Mirebeau on 29/08/13.
//
//

#ifndef ITKFM_RescaleAndExport_Function_h
#define ITKFM_RescaleAndExport_Function_h

#include "itkRescaleIntensityImageFilter.h"
#include "itkImageFileWriter.h"
#include "itkCastImageFilter.h"
#include <string>

// Simple function

template<typename InputImageType, typename OutputImageType = itk::Image<unsigned char, InputImageType::ImageDimension> >
bool RescaleAndExport(const InputImageType * input, std::string filename){
        
    typedef itk::ImageFileWriter<OutputImageType> WriterType;
    typename WriterType::Pointer writer = WriterType::New();

    typedef itk::RescaleIntensityImageFilter<InputImageType,OutputImageType> CasterType;
    typename CasterType::Pointer caster = CasterType::New();

    typedef itk::CastImageFilter<InputImageType, OutputImageType> CasterType2;
    typename CasterType2::Pointer caster2 = CasterType2::New();
    
    typedef typename OutputImageType::PixelType OutputPixelType;
    if(std::numeric_limits<OutputPixelType>::is_integer){
        caster->SetOutputMinimum(std::numeric_limits<OutputPixelType>::min());
        caster->SetOutputMaximum(std::numeric_limits<OutputPixelType>::max());
        caster->SetInput(input);
        writer->SetInput(caster->GetOutput());
    } else {
        caster2->SetInput(input);
        writer->SetInput(caster2->GetOutput());
    }
    
    writer->SetFileName(filename.c_str());
    
    try {
        writer->Update();
    } catch (itk::ExceptionObject & e) {
        std::cerr << "Error in Rescale and Export, when exporting image " << endl << input << endl << " to file " << filename << endl;
        std::cerr << e;
        return false;
    }
    
    return true;
}

namespace Tests {
    
    bool RescaleAndExport_Test(){
        typedef itk::Image<double,2> ImageType;
        
        typedef ImageType::RegionType RegionType;
        RegionType::SizeType size;
        size.Fill(10);
        RegionType::IndexType index;
        index.Fill(0);
        ImageType::RegionType region(index,size);
        
        auto image = ImageType::New();
        image->SetRegions(region);
        image->Allocate();
        for(index[0]=0; index[0]<(itk::IndexValueType)size[0]; ++index[0]){
            for(index[1]=0; index[1]<(itk::IndexValueType)size[1]; ++index[1]){
                image->GetPixel(index)
                = index[0]+index[1] < (itk::IndexValueType)size[0];
            } // for index[1]
        } // for index[0]
        
        return RescaleAndExport(image.GetPointer(), "RescaleAndExport_Test.bmp");
        
    } // void RescaleAndExport
    
} // anonymous namespace

#endif
