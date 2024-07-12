#ifndef tests_macro_hpp
    #include <string>
    #define tests_macro_hpp
    #ifndef VERBOSE
        #define VLOG(x) ((void) 0)
        #define VLOGe(x) ((void) 0)
    #endif 
    #ifdef VERBOSE
        #define VLOG(x) std::cout <<" [INFO] "<<x << std::endl
        #define VLOGe(x) std::cerr <<" [ERROR] "<<x << std::endl
    #endif 

    template <typename T>
    inline std::string toStringPers(const T& value) {
        return std::to_string(value);
    }
    template <>
    inline std::string toStringPers(const std::string& value) {
        return value;
    }
    #define CONV(x) toStringPers(x)
    #define CLOCK std::chrono::high_resolution_clock::now();
    #define COUNT_T(start,stop) std::chrono::duration_cast<std::chrono::milliseconds>(stop - start).count()
    #endif