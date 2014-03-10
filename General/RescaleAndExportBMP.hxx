//
//  RescaleAndExportBMP.hxx
//  itkDiffusion
//
//  Created by MIREBEAU Jean-Marie on 20/11/12.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#ifndef itkDiffusion_RescaleAndExportBMP_hxx 
#define itkDiffusion_RescaleAndExportBMP_hxx

namespace itk {
    
    template<typename TImage>
    RescaleAndExportBMP<TImage>
    ::RescaleAndExportBMP()
    {
        m_caster = CastFilterType::New();
        m_writer = WriterType::New();
        
        m_caster->SetOutputMinimum(0);
        m_caster->SetOutputMaximum(255);
        
        m_writer->SetInput(m_caster->GetOutput()); 
        m_filename = "default.bmp";
    }

    template<typename TImage>
    void
    RescaleAndExportBMP<TImage>
    ::Update()
    {
        m_writer->SetFileName(m_filename.c_str());
        m_writer->Update();
    }
    
}
#endif
