//
//  SyntheticImageSource.h
//  itkDiffusion
//
//  Created by Jean-Marie Mirebeau on 07/03/2014.
//
//

#ifndef itkDiffusion_SyntheticImageSource_h
#define itkDiffusion_SyntheticImageSource_h

#include "Macro.h"
#include "itkImageSource.h"
#include "itkImageRegionIteratorWithIndex.h"

namespace itk {
    template <typename TImage, typename TFunctor>
    class SyntheticImageSource : public ImageSource<TImage> {
    public:
        typedef SyntheticImageSource Self;
        typedef ImageSource<TImage> Superclass;
        typedef SmartPointer<Self> Pointer;
        typedef SmartPointer<const Self> ConstPointer;
        
        /// Method for creation through the object factory.
        itkNewMacro(Self);
        /// Run-time type information (and related methods).
        itkTypeMacro(SyntheticImageSource, Superclass);
        
        typedef TFunctor FunctorType;
        typedef TImage ImageType;
        static const int Dimension = ImageType::ImageDimension;
        typedef ImageBase<Dimension> PhysicalDataType;
        
        typedef double ScalarType;
        typedef Point<ScalarType,Dimension> PointType;
        
        GetSetFunctorMacro(Functor, FunctorType);
//        void SetFunctor(const FunctorType & f){m_Functor=f; this->Modified();}
//        itkGetMacro(Functor, FunctorType);
        
        itkGetModifiableObjectMacro(PhysicalData,PhysicalDataType);
//        itkGetMacro(PhysicalData, PhysicalDataType*);
        SyntheticImageSource(){m_PhysicalData = PhysicalDataType::New();}
        
    protected:
        FunctorType m_Functor;
        typename PhysicalDataType::Pointer m_PhysicalData;
        
        typedef typename ImageType::RegionType    OutputRegionType;
        
        virtual void GenerateOutputInformation(){
            if(this->GetOutput())
                this->GetOutput()->CopyInformation(m_PhysicalData);
        }
        
        virtual void ThreadedGenerateData(const OutputRegionType & region, ThreadIdType threadId){
            if(region.GetSize()[0]==0) return;
            
            ImageRegionIteratorWithIndex<ImageType> outputIt(this->GetOutput(), region);

            for(outputIt.GoToBegin(); !outputIt.IsAtEnd(); ++outputIt){
                PointType p;
                this->GetOutput()->TransformIndexToPhysicalPoint(outputIt.GetIndex(), p);
                outputIt.Set(m_Functor(p));
            }
        }
    };
}

#endif
