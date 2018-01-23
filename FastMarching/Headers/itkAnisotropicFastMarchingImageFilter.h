//
//  itkAFM.h
//  ITKFM
//
//  Created by MIREBEAU Jean-Marie on 08/11/12.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#ifndef ITKFM_itkAFM_h
#define ITKFM_itkAFM_h

#include "itkFastMarchingImageFilter_Modified.h"
#include "itkPolyLineParametricPath.h"
namespace itk {
    
    /** \class AnisotropicFastMarchingImageFilter
      \brief Solve an Anisotropic Eikonal equation using Anisotropic Fast Marching
     
     Anisotropic Fast Marching Image Filter
     
     Anisotropic Fast Marching solves an Eikonal equation where the speed
     is always non-negative and depends on the position and \a orientation.
     The speed can be anisotropic : e.g. up to 100x larger in one direction 
     than in another, without performance impact, and these directions need 
     not be aligned with the grid.
          
     Updates are preformed using an entropy satisfy scheme where only
     "upwind" neighborhoods are used. These neighborhoods are built using arithmetic 
     techniques, which guarantee a good accuracy and a low complexity of the filter, 
     almost independently of the speed anisotropy. 
     
     The construction and storage of these neighborhoods results in a CPU time 
     and memory usage roughly double that of the base class FastMarchingImageFilter
     (which is restricted to isotropic speeds).
     If speed is isotropic in your application, consider using the base class 
     FastMarchingImageFilter instead.

     Fast Marching sweeps through N grid points in \f${\mathcal O}(k N \log N)\f$ steps to obtain
     the arrival time value as the front propagates through the grid. The factor
     k depends logarithmically on the speed anisotropy.
     
     Implementation of this class is based on the papers
     Efficient Fast Marching with Finsler Metrics, J.-M. Mirebeau, 2012 (in 2D)
     Anisotropic Fast Marching using Lattice Basis Reduction, J.-M. Mirebeau, 2012 (in 3D)
     
     The input of the filter, provided through the method SetInput(), is a 
     metric : an image in which each point defines a norm, potentially anisotropic and asymmetric.
     This class is templated over the Norm type. Three types of norms are implemented so far:
     Riemannian2DNorm, Riemannian3DNorm and Finsler2DNorm. They allow to define Riemannian metrics in
     dimension 2, 3, and a special type of asymmetric metric in dimension 2.
    
     If the GenerateGradientImage flag is set, then the upwind gradient of image is saved, 
     in a raw form. The standard form of this gradient, an image of covariant vectors,
     can be obtained via the method GetGradientImage().

     If the GenerateGradientImage 
     flag is set, the shortest path from an arbitrary point of  
     to the closest seed can be extracted by the method Geodesic(). Geodesic extraction is 
     performed as in the second reference paper above.
     
     
     This class includes an additionnal template parameter internal_size_t, defaulted to unsigned int.
     It should not be modified unless the user plans to use images of resolution higher than 100 megapixels.
     In that case, change this parameter to size_t (this will increase memory usage). 
     
     Other relevant information is available in the documentation of the base class FastMarchingImageFilter
     
     * \ingroup LevelSetSegmentation
     * \ingroup ITKFastMarching
     */
    
    template <typename TNorm, typename internal_size_t = unsigned int>
    class ITK_EXPORT AnisotropicFastMarchingImageFilter : 
    public
    FastMarchingImageFilter_Modified<
    Image<typename TNorm::ValueType,TNorm::Dimension>,
    Image<TNorm,TNorm::Dimension>
    > {
    public:
        
        /** Standard class typdedefs. */
        typedef AnisotropicFastMarchingImageFilter    Self;
        typedef SmartPointer< Self >       Pointer;
        typedef SmartPointer< const Self > ConstPointer;
        
        typedef TNorm NormType;
        typedef typename NormType::ValueType ValueType;
        static const int Dimension = NormType::Dimension;
        
        typedef Image<ValueType,Dimension> TLevelSet;
        
        /** This class uses an image of norms, and not of speeds contrary 
         to the base class FastMarchingImageFilter. A speed value c corresponds to 
         the isotropic norm Identity / (c*c).
         */ 
        typedef Image<NormType,Dimension> NormImageType;
        typedef FastMarchingImageFilter_Modified<TLevelSet,NormImageType> Superclass;
        
        
        typedef Index<Dimension> IndexType;
        typedef typename NormType::VectorType VectorType;
        typedef itk::Point<ValueType,Dimension> PointType;
        
        typedef ImageRegion<Dimension> ImageRegionType;
        typedef typename Superclass::LevelSetImageType LevelSetImageType;

        /// Method for creation through the object factory. 
        itkNewMacro(Self);
        
        /// Run-time type information (and related methods). 
        itkTypeMacro(AnisotropicFastMarchingImageFilter, FastMarchingImageFilter);
                
        /** The GenerateGradientImage flag must be set for the methods
         GetGradientImage() AND Geodesic() to work.
         */
        virtual bool GetGenerateUpwindGradient() const {return this->m_GenerateUpwindGradient;}
        virtual void SetGenerateUpwindGradient(bool Generate) {this->m_GenerateUpwindGradient = Generate;}

        typedef CovariantVector<ValueType,Dimension> CovariantVectorType;
        typedef Image<CovariantVectorType,Dimension> GradientImageType;
        /// Returns the upwind gradient of the output.
        typename GradientImageType::Pointer ComputeGradientImage();
        
        /** The input vector to the Geodesic method should contain a single point,
         given by its GeodesicContinuousIndex : the starting point of the minimal path that 
         the user wants to extract.
         After execution, this vector contains a sequence of GeodesicContinuousIndices (p_i,v_i). The
         values of the distance at indices p_i are strictly decreasing, and the last one is a 
         seed formerly provided by the user.
         */
        
        typedef ContinuousIndex<ValueType,Dimension> ContinuousIndexType;
        
        bool Geodesic(std::vector<ContinuousIndexType> &) const;
        bool Geodesic(std::vector<PointType> &) const;
        
    protected:
        
        /** A GeodesicContinuousIndex is a pair containing a proper index,
         and a (small) additive correction allowing to point off the grid. Conversions to and from
         itk::ContinuousIndex are provided.
         */
        
        struct GeodesicContinuousIndex : public std::pair<IndexType,VectorType> {
            ContinuousIndexType GetContinuousIndex() const;
            template<typename TComponent> GeodesicContinuousIndex(const ContinuousIndex<TComponent,Dimension> &);
            GeodesicContinuousIndex(const IndexType &p, const VectorType &v):std::pair<IndexType,VectorType>(p,v){};
        };
        
        typedef std::vector<GeodesicContinuousIndex> PathRawType;
        void Geodesic(PathRawType &) const;
        
        
        AnisotropicFastMarchingImageFilter(){};
        ~AnisotropicFastMarchingImageFilter(){};
        
        virtual void Initialize(LevelSetImageType *) override;
        
        virtual void UpdateNeighbors(const IndexType & index,
                                     const NormImageType *, LevelSetImageType *) override;
        
        /** Due to the inherited structure, I have to pass the two following arguments from UpdateNeighbors to 
         UpdateValue as class range variables */
        IndexType YoungestAliveIndex;
        unsigned int YoungestAliveIndex_i;
        
        virtual double UpdateValue(const IndexType & index,
                                   const NormImageType *, LevelSetImageType *) override;
        
        bool m_GenerateUpwindGradient;
        
        struct GradientRawData {
            typedef typename NormType::IndicesType IndicesType;
            typedef typename NormType::WeightsType WeightsType;
            
            typedef std::pair<IndicesType, WeightsType> RawType;
            typedef Image<RawType,Dimension> RawImageType;
            typename RawImageType::Pointer Raw;
            
            void Init(ImageRegionType);        
        } Gradient; 

        
        typedef typename NormType::ShortOffsetType ShortOffsetType;
        typedef typename NormType::Stencil DirectStencilType;
        
        struct StencilData {
            
            ImageRegionType BufferedRegion;

            typedef typename ShortOffsetType::ValueType ShortOffsetValueType;
            typedef typename NormType::CompressedOffsetType CompressedOffsetType;

            
            /// The start and size of the CompressedStencil associated to a pixel, given by its index.
            DirectStencilType Direct(const IndexType &P) const;
            
            typedef internal_size_t StencilStartType;
            
            class ReverseStencilType {                
                const StencilStartType size_;
                const ShortOffsetType * const start_;
                const ShortOffsetValueType * const start_i_;
            public:
                StencilStartType size() const {return size_;}
                const ShortOffsetType & operator[](int i) const {return start_[i];}
                ShortOffsetValueType _i(int i) const {return start_i_[i];}
                ReverseStencilType(StencilStartType size, ShortOffsetType * start, ShortOffsetValueType * start_i):size_(size), start_(start), start_i_(start_i){};
            };
                        
            /// The start position and size of the ReversedStencil associated to a pixel, given by its index.
            const ReverseStencilType Reverse(const IndexType &P);
            
            /// Initialization of the stencils. Calls the private methods Allocate, SetCompressed and Set Reversed.
            void Init(const NormImageType*);
            
        protected:
            void Allocate();
            void SetCompressed(const NormImageType*); 
            void SetReversed();
            
            
            /** The compressed stencils, associated to the norms in the input image, are arranged linearly 
             in a large array. They may have different sizes. The CompressedStart image keeps track of 
             of the start of each stencil. Compressed stencils are "uncompressed", yielding the direct stencils,
             by calling the GetStencil method of the NormType. */
            std::vector<CompressedOffsetType> Compressed;
            
            /** The reversed stencils. u is in the reversed stencil of v, if and only if v is in the direct 
             stencil of u. Arranged linearly. The reversed start image keeps track of the start. */
            std::vector<ShortOffsetType> Reversed;
            
            /** If v is the i-th element in the stencil of u, then u is in the reversed stencil of v, and the 
             integer i is kept in the following array. */
            std::vector<ShortOffsetValueType> Reversed_i;
            
            /** The positions of the consecutive stencils, in the Compressed and the Reversed arrays, are kept 
             in the CompressedStart and the ReverseStart images. */
            typedef Image<StencilStartType,Dimension> StartImageType;
            typename StartImageType::Pointer CompressedStart;
            typename StartImageType::Pointer ReversedStart;

#pragma remark("Eliminate these two iterators ?")
            mutable ImageRegionIterator<StartImageType> CompressedStartIt;
            ImageRegionIterator<StartImageType> ReversedStartIt;
            friend typename GradientImageType::Pointer AnisotropicFastMarchingImageFilter::GetGradientImage();
        } Stencils;
    };
    
}

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkAnisotropicFastMarchingImageFilter.hxx"
#include "itkAnisotropicFastMarchingImageFilter_StencilData.hxx"
#include "itkAnisotropicFastMarchingImageFilter_Geodesic.hxx"
#endif

#endif
