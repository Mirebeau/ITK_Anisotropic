//
//  AnisotropicDiffusionLBRImageFilter.hxx
//  itkDiffusion
//
//  Created by Jean-Marie Mirebeau on 28/02/2014.
//
//

#ifndef itkDiffusion_AnisotropicDiffusionLBRImageFilter_hxx
#define itkDiffusion_AnisotropicDiffusionLBRImageFilter_hxx

namespace itk {

    template<typename TI>
    void AnisotropicDiffusionLBRImageFilter<TI>::GenerateData(){
        
        {
            typedef CastImageFilter<ImageType, ScalarImageType> CasterType;
            typename CasterType::Pointer caster = CasterType::New();
            caster->SetInput(this->GetInput());
            caster->Update();
            image = caster->GetOutput();
        }
        
        while(m_DiffusionTime>0){
            ComputeDiffusionTensors();
            linearDiffusionFilter->SetInputImage(image);
            linearDiffusionFilter->SetInputTensor(tensorImage);
            linearDiffusionFilter->SetMaxDiffusionTime(m_DiffusionTime);
            linearDiffusionFilter->Update();
            image = linearDiffusionFilter->GetOutput();
            m_DiffusionTime -= linearDiffusionFilter->GetEffectiveDiffusionTime();
            
            m_LinearFilterEffectiveTimesAndIterations.push_back(std::pair<ScalarType,int>(linearDiffusionFilter->GetEffectiveDiffusionTime(),linearDiffusionFilter->GetEffectiveNumberOfTimeSteps()));
        }
        
        {
            typedef CastImageFilter<ScalarImageType, ImageType> CasterType;
            typename CasterType::Pointer caster = CasterType::New();
            caster->SetInput(image);
            caster->Update();
            this->GraftOutput(caster->GetOutput());
        }
    }
    
    template<typename TI>
    void AnisotropicDiffusionLBRImageFilter<TI>::ComputeDiffusionTensors(){
        structureTensorFilter->SetInput(image);
        
        struct DiffusionTensorFunctor {
            Self * eigenValuesFunctor;
            TensorType operator()(const TensorType & S){
                EigenValuesArrayType eigenValues;
                typename TensorType::EigenVectorsMatrixType eigenVectors;
                S.ComputeEigenAnalysis(eigenValues,eigenVectors);
                
                // For convenience, eigenvalues are sorted
                Vector<int,Dimension> order;
                for(int i=0; i<Dimension; ++i) order[i]=i;
                
                struct OrderingType { // c++ 11 would be : [& eigenValues](int i, int j)->bool {return eigenValues[i]<eigenValues[j];}
                    bool operator()(int i, int j) const {return this->e[i]<this->e[j];}
                    const EigenValuesArrayType & e;
                    OrderingType(const EigenValuesArrayType & e_):e(e_){};
                } ordering(eigenValues);
                
                std::sort(order.Begin(), order.End(),ordering);
                
                
                std::sort(eigenValues.Begin(), eigenValues.End());
                EigenValuesArrayType ev = this->eigenValuesFunctor->EigenValuesTransform(eigenValues);
                
                TensorType DiffusionTensor;
                for(int i=0; i<Dimension; ++i){
                    DiffusionTensor(order[i],order[i]) = ev[i];
                    for(int j=0; j<i; ++j) DiffusionTensor(i,j) = 0.;
                }
                
                return DiffusionTensor.Rotate(eigenVectors.GetTranspose());
            }
        };
        
        typedef UnaryFunctorImageFilter<TensorImageType, TensorImageType, DiffusionTensorFunctor> ImageFunctorType;
        auto imageFunctor = ImageFunctorType::New();
        imageFunctor->GetFunctor().eigenValuesFunctor = this;
        imageFunctor->SetInput(structureTensorFilter->GetOutput());
        
        imageFunctor->Update();
        tensorImage=imageFunctor->GetOutput();        
    }
    
}




#endif
