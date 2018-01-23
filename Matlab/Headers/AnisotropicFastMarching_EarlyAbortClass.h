//
//  AnisotropicFastMarching_Vijaya.h
//  MatlabITK
//
//  Created by Jean-Marie Mirebeau on 16/09/13.
//
//

#ifndef MatlabITK_AnisotropicFastMarching_Vijaya_h
#define MatlabITK_AnisotropicFastMarching_Vijaya_h

#include "itkAnisotropicFastMarchingImageFilter.h"
#include "mexMessageWrapper.h"

namespace itk {
    
    template<typename TNorm>
    struct AnisotropicFastMarchingImageFilter_EarlyAbort :
    public AnisotropicFastMarchingImageFilter<TNorm> {
        
        // ITK Stuff
        typedef AnisotropicFastMarchingImageFilter_EarlyAbort    Self;
        typedef AnisotropicFastMarchingImageFilter<TNorm> SuperClass;

        typedef SmartPointer< Self >       Pointer;
        typedef SmartPointer< const Self > ConstPointer;
        itkNewMacro(Self);
        itkTypeMacro(AnisotropicFastMarchingImageFilter_EarlyAbort, AnisotropicFastMarchingImageFilter);

        // required typedefs
        static const unsigned int Dimension = SuperClass::Dimension;
        typedef Index<Dimension> IndexType;
        typedef typename SuperClass::NormType NormType;
        typedef typename SuperClass::ValueType ValueType;
        typedef typename SuperClass::LevelSetImageType LevelSetImageType;
        typedef typename SuperClass::VectorType VectorType;
        typedef typename SuperClass::PointType PointType;
        typedef typename SuperClass::Superclass::OutputSpacingType SpacingType;
        
        void StopWhenLastIsAccepted(const IndexType & indexToAccept){
            selectedStoppingCriteria[LastAccepted] = true;
            m_RemainingIndicesToAccept.insert(indexToAccept);
        }
        
        void StopWhenFirstIsAccepted(const IndexType & indexToAccept){
            selectedStoppingCriteria[FirstAccepted] = true;
            m_StopWhenFirstIsAccepted.insert(indexToAccept);
        }
        
        void StopAtDistance(ValueType distance){
            selectedStoppingCriteria[Distance] = true;
            m_Distance = distance;
        }
        
        void StopAtEuclideanDistance(ValueType euclideanDistance){
            selectedStoppingCriteria[EuclideanDistance] = true;
            m_EuclideanDistance = euclideanDistance;
        }
        
        typename LevelSetImageType::Pointer GetEuclideanPathLengths() const {return EuclideanPathLengths;}
        void SetEuclideanSpacing(SpacingType s){m_EuclideanSpacing=s;}
        
        enum StoppingCriteria {LastAccepted=0, FirstAccepted, Distance, EuclideanDistance, None, nStoppingCriteria};
        
        StoppingCriteria GetActiveStoppingCriteria() const {return activeStoppingCriteria;}
        
        IndexType GetStoppingIndex() const {return stoppingIndex;}

    protected:
        bool selectedStoppingCriteria[nStoppingCriteria] = {false,false,false,false,false};
        StoppingCriteria activeStoppingCriteria;
        IndexType stoppingIndex;
        
        
        std::set<IndexType, typename IndexType::LexicographicCompare> m_RemainingIndicesToAccept, m_StopWhenFirstIsAccepted;
        
        ValueType m_Distance;
        
        ValueType m_EuclideanDistance;
        typename LevelSetImageType::Pointer EuclideanPathLengths;
        SpacingType m_EuclideanSpacing;
        
        typedef typename SuperClass::NormImageType NormImageType;
        
        virtual void Initialize(LevelSetImageType * output) override {
            activeStoppingCriteria = None;
            SuperClass::Initialize(output);
            
            if(selectedStoppingCriteria[EuclideanDistance]){
                if(!this->GetGenerateUpwindGradient())
                    itkExceptionMacro("StopAtEuclideanDistance needs the gradients");
                EuclideanPathLengths=LevelSetImageType::New();
                EuclideanPathLengths->SetRegions(this->GetOutputRegion());
                EuclideanPathLengths->CopyInformation(output);
                EuclideanPathLengths->Allocate();
                EuclideanPathLengths->FillBuffer(std::numeric_limits<ValueType>::max());
            }
        }

        
        // Called when point "index" is accepted (i.e. frozen)
        virtual void UpdateNeighbors(const IndexType & index,
                                     const NormImageType * norms,
                                     LevelSetImageType * distance) override {
            if(activeStoppingCriteria != None) return;
            stoppingIndex = index;
            
            // Vijaya's stopping criteria
            if(selectedStoppingCriteria[Distance]){
                if(distance->GetPixel(index) > m_Distance){
                    activeStoppingCriteria = Distance;
                    return;
                }
            }
            
            if(selectedStoppingCriteria[LastAccepted]){
                if(m_RemainingIndicesToAccept.empty()){
                    activeStoppingCriteria = LastAccepted;
                    return;
                }
                m_RemainingIndicesToAccept.erase(index);
            }
            
            if(selectedStoppingCriteria[FirstAccepted]){
                if(m_StopWhenFirstIsAccepted.find(index) != m_StopWhenFirstIsAccepted.end()){
                    activeStoppingCriteria = FirstAccepted;
                    return;
                }
            }
            
            // Chenda's stopping criteria
            if(selectedStoppingCriteria[EuclideanDistance]){
                const auto & grad = this->Gradient.Raw->GetPixel(index);
                const auto & indices = grad.first;
                const auto & weights = grad.second;
                const auto Stencil = this->Stencils.Direct(index);
                
                auto & euclideanDistance = EuclideanPathLengths->GetPixel(index);
                
                euclideanDistance=0.;
                if(!weights.IsNullWeights()) {
                    VectorType u(0.);
                    
                    for(int k=0; k<Dimension; ++k){
                        if(indices[k] == NormType::IndexToIgnore())
                            continue;
                        
                        const ValueType w = weights(k);
                        const auto & offsetNeigh = Stencil[indices[k]];
                        IndexType indexNeigh;
                        for(int l=0; l<Dimension; ++l)
                            indexNeigh[l] = index[l]+offsetNeigh[l];
                        euclideanDistance += w*EuclideanPathLengths->GetPixel(indexNeigh);
                        
                        for(int l=0; l<Dimension; ++l)
                            u[l] += w*offsetNeigh[l]*m_EuclideanSpacing[l];
                        
                    }
                    euclideanDistance += u.GetNorm();
                }
                if(euclideanDistance >= m_EuclideanDistance){
                    activeStoppingCriteria = EuclideanDistance;
                    return;
                }
            }
            
            // Perform as usual
            SuperClass::UpdateNeighbors(index, norms, distance);
        }
    };

    
}

#endif
