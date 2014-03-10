//
//  UnrollingConstructors.h
//  MatlabITKFastMarching
//
//  Created by Jean-Marie Mirebeau on 05/04/13.
//
//

#ifndef MatlabITKFastMarching_UnrollingConstructors_h
#define MatlabITKFastMarching_UnrollingConstructors_h


namespace UnrollingConstructors {
    
    template<int i, typename ContainerType>
    void Set(ContainerType & V){
        static_assert(i==ContainerType::Dimension,"Wrong number of arguments");
        return;
    }
    
    template<int i=0, typename ContainerType, typename T, typename... Args>
    void Set(ContainerType & V, T x, const Args&... xx){
        V[i]=x;
        //(*V) = x;
        Set<i+1>(V, xx...);
    }
    
    
    template<typename ContainerType, typename... Args>
    ContainerType Make(const Args&... xx){
        ContainerType Container;
        Set(Container,xx...);
        return Container;
    }
    
    /*
    
    void printf(const char* s) {
        while (*s) {
            if (*s == '%' && *++s != '%')
                throw std::runtime_error("invalid format string: missing arguments");
            std::cout << *s++;
        }
    }
    
    template<typename T, typename... Args>
    void printf(const char* s, const T& value, const Args&... args) {
        while (*s) {
            if (*s == '%' && *++s != '%') {
                std::cout << value;
                return printf(++s, args...);
            }
            std::cout << *s++;
        }
        throw std::runtime_error("extra arguments provided to printf");
    }
    */
    
    
    /*
    namespace __UnrollingConstructors__ {
        
        
        template<typename ContainerType, typename ValueType, int i, int Dimension>
        void Set(ContainerType & Container, ValueType Value){
//            static_assert(i==Dimension-1,"Unrolling constructor error : too few arguments provided.");
            Container[i] = Value;
        }

        
        template<typename ContainerType, typename T, typename ...ValueType, int i, int Dimension>
        void Set(ContainerType & Container, T FirstValue, ValueType ...RemainingValues){
            if(i==Dimension) return;
            static_assert(i<Dimension,"Unrolling constructor error : too many arguments provided.");
            Container[i] = FirstValue;
            Set<ContainerType, T, ValueType..., i+1, Dimension>(Container, RemainingValues...);
        }
    }
    

    template<typename ContainerType, typename ValueType = typename ContainerType::ValueType, typename ValuesType = ...ValueType, int Dimension = ContainerType::Dimension>
    ContainerType Make(ValueTypes ...Values){
        ContainerType Container;
        __UnrollingConstructors__::Set<ContainerType,ValueType,0,Dimension>(Container,Values...);
        return Container;
    }
 
    
    template<typename ContainerType, typename T, typename ...Args, int Dimension = ContainerType::Dimension>
    ContainerType Make(T Value, Args ...Values){
        ContainerType Container;
        __UnrollingConstructors__::Set<ContainerType, T, Args...,0,Dimension>(Container,Value, Values...);
        return Container;
    }
    */
    
}


#endif
