//
//  ExportToMathematica.h
//  ITKFM
//
//  Created by Jean-Marie Mirebeau on 10/01/13.
//
//

#ifndef ITKFM_ExportToMathematica_h
#define ITKFM_ExportToMathematica_h

#include <sstream>
#pragma remark("Get rid of this file, use HDF5 and Print Range As List instead")
namespace itk {
    
    template<typename ForwardIteratorType>
    void PrintRangeAsList(std::ostream & f, ForwardIteratorType begin, ForwardIteratorType end){
        f << "{";
        if(begin==end) {f << "}"; return;}
        for(f << *begin++; begin!=end; ++begin)
            f << "," << *begin;
        f << "}";
    }
    
    
    // ******** For now, 3D images, and geodesics, are exported in Mathematica friendly format as follows *********
    template<class ImageIterator, unsigned int Dimension>
    void ExportTensorToTXT(std::ostream &f, ImageIterator it, Size<Dimension> N, bool one_per_line=false)
    {
        it.GoToBegin();
        const unsigned int dim = N.Dimension;
        int prod[dim]; prod[0]=N[0]; for(int k=1; k<dim; ++k) prod[k]=N[k]*prod[k-1];
        for(int k=0; k<dim; ++k) f<<"{";
        if(!it.IsAtEnd()) {f<<it.Value(); ++it;}
        int i=1;
        for(; !it.IsAtEnd(); ++it, ++i){
            for(int k=0; k<dim; ++k) if(i%(prod[k])==0) f<<"}";
            f<<","; if(one_per_line) f<<"\n";
            for(int k=0; k<dim; ++k) if(i%(prod[k])==0) f<<"{";
            f << it.Value();
        }
        for(int k=0; k<dim; ++k) f<<"}";
        assert(i==prod[dim-1]); //size and format must match
    }

    
#pragma message("To do : should not be used anymore, thanks to the new geodesic functions")
    template<class GeodesicContinuousIndex>
    void ExportGeodesicToTXT(std::ostream &f, std::vector<GeodesicContinuousIndex> &geodesic)
    {
        const unsigned int Dimension = GeodesicContinuousIndex::first_type::Dimension;
        f<<"{";
        for(int i=0; i<geodesic.size(); ++i){
            if(i>0) f << "," << endl;
            f << "{";
            for(int j=0; j<Dimension; ++j) {
                if(j>0) f << ",";
                f << geodesic[i].first[j]+geodesic[i].second[j];
            }
            f << "}";
        }
        f << "}";
    }
    
    
    class MathematicaExporterType {
        std::ofstream & f; // Cannot be a data member due to implicit copies in std::for_each
        bool isFirstWrite=true;
            
        // Special format for FixedArray, for compliance with Mathematica
        template<typename T, unsigned int D>
        void Print(const itk::Point<T,D> & p){Print(static_cast<FixedArray<T,D> >(p));}
        template<typename T, unsigned int D>
        void Print(const itk::ContinuousIndex<T,D> & p){Print(static_cast<FixedArray<T,D> >(p));}
        
        template<typename T, unsigned int D>
        void Print(const itk::FixedArray<T,D> & p){
            assert(D>0); // Only supports arrays of positive length
//            static_assert(D>0,"Only supports arrays of positive length");
            f << "{" << p[0];
            for(int i=1; i<D; ++i)
                f << "," << p[i];
            f << "}";
        }
        
        template<typename DataType>
        void Print(const DataType & p){f << p;}
        
    public:
        void close(){if(f.is_open()) {f << "}"; f.close();} }
        void open(std::string filename){
            close();
            f.open(filename.c_str());
            f << "{";
            isFirstWrite=true;
        }
        
        MathematicaExporterType(std::ofstream & f_):f(f_){};
        ~MathematicaExporterType(){close();}

        template<typename DataType>
        void operator()(const DataType & p)
        {if(isFirstWrite) isFirstWrite=false; else f << ","; Print(p);}
    };
            
    
} // namespace itk

#endif
