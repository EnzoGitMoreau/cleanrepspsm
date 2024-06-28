// Cytosim was created by Francois Nedelec. Copyright 2020 Cambridge University.
// FJN, Strasbourg 12.06.2018

#ifndef MATRIX44
#define MATRIX44

#include "real.h"
#include "vector3.h"
#include "vector4.h"
#include <cstdio>
#include <iostream>

#if defined(__AVX__) && REAL_IS_DOUBLE
#  include "simd.h"
#  define MATRIX44_USES_AVX 1
#elif defined(__SSE3__) && !REAL_IS_DOUBLE
#  include "simd_float.h"
#  define MATRIX44_USES_AVX 0
#else
#  define MATRIX44_USES_AVX 0
#endif


/// 4x4 matrix class with 16 'real' elements stored in column order
class alignas(4*sizeof(real)) Matrix44
{
public:
    
    /// unsigned integer type used for indices
    typedef size_t index;

    /// access to modifiable element by index
    real& operator[](index i)       { return val[i]; }
    
    /// access element value by index
    real  operator[](index i) const { return val[i]; }

public:

    /// values of the elements
    real val[4*4];
    
    Matrix44() {}
    
    /// copy constructor
    Matrix44(Matrix44 const& M)
    {
        for ( index u = 0; u < 16; ++u )
            val[u] = M.val[u];
    }
    
    /// construct Matrix from coordinates (column-major)
    Matrix44(real a, real b, real c, real d,
             real e, real f, real g, real h,
             real i, real j, real k, real l,
             real m, real n, real o, real p )
    {
        val[0x0] = a;
        val[0x1] = b;
        val[0x2] = c;
        val[0x3] = d;
        val[0x4] = e;
        val[0x5] = f;
        val[0x6] = g;
        val[0x7] = h;
        val[0x8] = i;
        val[0x9] = j;
        val[0xA] = k;
        val[0xB] = l;
        val[0xC] = m;
        val[0xD] = n;
        val[0xE] = o;
        val[0xF] = p;
    }

    /// construct Matrix with diagonal terms set to `d` and other terms set to `z`
    Matrix44(real z, real d)
    {
        val[0x0] = d;
        val[0x1] = z;
        val[0x2] = z;
        val[0x3] = z;
        val[0x4] = z;
        val[0x5] = d;
        val[0x6] = z;
        val[0x7] = z;
        val[0x8] = z;
        val[0x9] = z;
        val[0xA] = d;
        val[0xB] = z;
        val[0xC] = z;
        val[0xD] = z;
        val[0xE] = z;
        val[0xF] = d;
    }

    ~Matrix44() {}
    
#pragma mark -
        
    /// dimensionality
    static constexpr size_t dimension() { return 4; }
    
    /// human-readable identifier
    static std::string what() { return "16"; }

    /// set all elements to zero
    void reset()
    {
        for ( index u = 0; u < 16; ++u )
            val[u] = 0.0;
    }
    
    /// set diagonal to 'dia' and other elements to 'off'
    void reset1(real off, real dia)
    {
        for ( index u = 0; u < 16; ++u )
            val[u] = off;
        for ( index u = 0; u < 4; ++u )
            val[u*5] = dia;
    }

    /// true if any value is different from 'zero'
    bool operator != (real zero) const
    {
        for ( index u = 0; u < 16; ++u )
            if ( val[u] != zero )
                return true;
        return false;
    }
    
    /// conversion to pointer of real
    //operator real const*() const { return val; }

    /// return modifiable pointer of 'real'
    real* data() { return val; }
    
#if defined(__AVX__) && REAL_IS_DOUBLE
    vec4 data0() const { return streamload4(val); }
    vec4 data1() const { return streamload4(val+4); }
    vec4 data2() const { return streamload4(val+8); }
    vec4 data3() const { return streamload4(val+12); }
#elif defined(__SSE3__) && !REAL_IS_DOUBLE
    vec4f data0() const { return streamload4f(val); }
    vec4f data1() const { return streamload4f(val+4); }
    vec4f data2() const { return streamload4f(val+8); }
    vec4f data3() const { return streamload4f(val+12); }
#endif

    /// unmodifiable pointer of real
    real const* data() const { return val; }

    /// address of element at line i, column j
    real* addr(const index i, const index j) { return val + ( i + 4*j ); }
    /// value of element at line i, column j
    real value(const index i, const index j) const { return val[i+4*j]; }

    /// access functions to element by line and column indices
    real& operator()(const index i, const index j)       { return val[i+4*j]; }
    real  operator()(const index i, const index j) const { return val[i+4*j]; }
    
    /// set elements from given array
    void load(const real ptr[])
    {
        for ( index u = 0; u < 16; ++u )
            val[u] = ptr[u];
    }

    /// copy elements to given array
    void store(real ptr[]) const
    {
        for ( index u = 0; u < 16; ++u )
            ptr[u] = val[u];
    }

    /// extract column vector at given index
    Vector4 column(const index i) const
    {
        return Vector4(val+4*i);
    }
    
    /// extract line vector at given index
    Vector4 line(const index i) const
    {
        return Vector4(val[i], val[4+i], val[8+i], val[12+i]);
    }
    
    /// extract diagonal
    Vector4 diagonal() const
    {
        return Vector4(val[0], val[5], val[10], val[15]);
    }
    
    /// sum of diagonal terms
    real trace() const
    {
        return ( val[0x0] + val[0x5] + val[0xA] + val[0xF] );
    }

    /// output in human-friendly format
    void print(FILE * f) const
    {
        fprintf(f, " / %9.3f %+9.3f %+9.3f %+9.3f \\\n", val[0x0], val[0x4], val[0x8], val[0xC]);
        fprintf(f, " | %9.3f %+9.3f %+9.3f %+9.3f |\n" , val[0x1], val[0x5], val[0x9], val[0xD]);
        fprintf(f, " | %9.3f %+9.3f %+9.3f %+9.3f |\n" , val[0x2], val[0x6], val[0xA], val[0xE]);
        fprintf(f, " \\ %9.3f %+9.3f %+9.3f %+9.3f /\n", val[0x3], val[0x7], val[0xB], val[0xF]);
    }
    
    /// print [ line1; line2; line3 ]
    void print(std::ostream& os) const
    {
        const int w = (int)os.width();
        os << std::setw(1) << "[";
        for ( index i = 0; i < 4; ++i )
        {
            for ( index j = 0; j < 4; ++j )
                os << " " << std::fixed << std::setw(w) << value(i,j);
            if ( i < 3 )
                os << ";";
            else
                os << " ]";
        }
    }

    /// conversion to string
    std::string to_string(std::streamsize w, std::streamsize p) const
    {
        std::ostringstream os;
        os.precision(p);
        os.width(w);
        print(os);
        return os.str();
    }

    /// scale all elements
    void scale(const real alpha)
    {
        for ( index u = 0; u < 16; ++u )
            val[u] *= alpha;
    }
    
    /// scale matrix
    void operator *=(const real alpha)
    {
        scale(alpha);
    }
    
    /// return opposite matrix (i.e. -M)
    const Matrix44 operator -() const
    {
        Matrix44 M;
        for ( index u = 0; u < 16; ++u )
            M.val[u] = -val[u];
        return M;
    }
    
    /// returns alpha * M
    const Matrix44 operator *(const real alpha) const
    {
        Matrix44 M;
        for ( index u = 0; u < 16; ++u )
            M.val[u] = val[u] * alpha;
        return M;
    }
    
    /// multiplication by scalar
    friend const Matrix44 operator *(const real alpha, Matrix44 const& mat)
    {
        return mat * alpha;
    }

    /// return sum of two matrices
    const Matrix44 operator +(Matrix44 const& M) const
    {
        Matrix44 res;
        for ( index u = 0; u < 16; ++u )
            res.val[u] = val[u] + M.val[u];
        return res;
    }
    
    /// return difference of two matrices
    const Matrix44 operator -(Matrix44 const& M) const
    {
        Matrix44 res;
        for ( index u = 0; u < 16; ++u )
            res.val[u] = val[u] - M.val[u];
        return res;
    }

    /// subtract given matrix
    void operator +=(Matrix44 const& M)
    {
#if MATRIX44_USES_AVX
        store4(val   , add4(load4(val   ), load4(M.val   )));
        store4(val+4 , add4(load4(val+4 ), load4(M.val+4 )));
        store4(val+8 , add4(load4(val+8 ), load4(M.val+8 )));
        store4(val+12, add4(load4(val+12), load4(M.val+12)));
#else
        for ( index u = 0; u < 16; ++u )
            val[u] += M.val[u];
#endif
    }

    /// add given matrix
    void operator -=(Matrix44 const& M)
    {
#if MATRIX44_USES_AVX
        store4(val   , sub4(load4(val   ), load4(M.val   )));
        store4(val+4 , sub4(load4(val+4 ), load4(M.val+4 )));
        store4(val+8 , sub4(load4(val+8 ), load4(M.val+8 )));
        store4(val+12, sub4(load4(val+12), load4(M.val+12)));
#else
        for ( index u = 0; u < 16; ++u )
            val[u] -= M.val[u];
#endif
    }
    
    /// transpose matrix in place
    void transpose()
    {
        std::swap(val[0x4], val[0x1]);
        std::swap(val[0x8], val[0x2]);
        std::swap(val[0xC], val[0x3]);
        std::swap(val[0x9], val[0x6]);
        std::swap(val[0xD], val[0x7]);
        std::swap(val[0xE], val[0xB]);
    }
    
    /// return transposed matrix
    Matrix44 transposed() const
    {
        Matrix44 res;
#if MATRIX44_USES_AVX
        vec4 v0 = load4(val);
        vec4 v1 = load4(val+4);
        vec4 v2 = load4(val+8);
        vec4 v3 = load4(val+12);
        vec4 u0 = unpacklo4(v0, v1);
        vec4 u1 = unpackhi4(v0, v1);
        v0 = unpacklo4(v2, v3);
        v1 = unpackhi4(v2, v3);
        v2 = catshift2(u0, v0);
        v3 = catshift2(u1, v1);
        store4(res.val   , blend22(u0, v2));
        store4(res.val+4 , blend22(u1, v3));
        store4(res.val+8 , blend22(v2, v0));
        store4(res.val+12, blend22(v3, v1));
#else
        for ( index x = 0; x < 4; ++x )
        for ( index y = 0; y < 4; ++y )
            res.val[y+4*x] = val[x+4*y];
#endif
        return res;
    }
    
    /// return scaled transposed matrix
    Matrix44 transposed(real alpha) const
    {
        Matrix44 res;
#if MATRIX44_USES_AVX
        vec4 a = set4(alpha);
        vec4 v0 = mul4(a, load4(val));
        vec4 v1 = mul4(a, load4(val+4));
        vec4 v2 = mul4(a, load4(val+8));
        vec4 v3 = mul4(a, load4(val+12));
        vec4 u0 = unpacklo4(v0, v1);
        vec4 u1 = unpackhi4(v0, v1);
        v0 = unpacklo4(v2, v3);
        v1 = unpackhi4(v2, v3);
        v2 = catshift2(u0, v0);
        v3 = catshift2(u1, v1);
        store4(res.val   , blend22(u0, v2));
        store4(res.val+4 , blend22(u1, v3));
        store4(res.val+8 , blend22(v2, v0));
        store4(res.val+12, blend22(v3, v1));
#else
        for ( index x = 0; x < 4; ++x )
        for ( index y = 0; y < 4; ++y )
            res.val[y+4*x] = alpha * val[x+4*y];
#endif
        return res;
    }

    /// maximum of all component's absolute values
    real norm_inf() const
    {
        real res = abs_real(val[0]);
        for ( index i = 1; i < 16; ++i )
            res = std::max(res, abs_real(val[i]));
        return res;
    }

    /// copy values from lower triangle to upper triangle
    void copy_lower()
    {
        val[0x4] = val[0x1];
        val[0x8] = val[0x2];
        val[0xC] = val[0x3];
        val[0x9] = val[0x6];
        val[0xD] = val[0x7];
        val[0xE] = val[0xB];
    }

    /// relative asymmetry of matrix (divided by the trace)
    real asymmetry() const
    {
        real t = abs_real(val[0]) + abs_real(val[0x5]) + abs_real(val[0xA]) + abs_real(val[0xF]);
        return (  abs_real(val[0x4]-val[0x1])
                + abs_real(val[0x8]-val[0x2])
                + abs_real(val[0xC]-val[0x3])
                + abs_real(val[0x9]-val[0x6])
                + abs_real(val[0xD]-val[0x7])
                + abs_real(val[0xE]-val[0xB]) ) / t;
    }
    
#pragma mark -

#if MATRIX44_USES_AVX
    
    /// multiplication by a 3-components vector: this * V = { x, y, z, garbage }
    const vec4 vecmul3_avx(const vec4 xyzt) const
    {
        vec4 xyxy = duplo2f128(xyzt);
        vec4 ztzt = duphi2f128(xyzt);
        vec4 xxxx = duplo4(xyxy);
        vec4 yyyy = duphi4(xyxy);
        vec4 zzzz = duplo4(ztzt);
        xxxx = mul4(load4(val), xxxx);
        yyyy = mul4(load4(val+4), yyyy);
        return fmadd4(load4(val+8), zzzz, add4(xxxx, yyyy));
    }

    /// multiplication by a 3-components vector: transpose(M) * V
    const vec4 trans_vecmul3_avx(double const* V) const
    {
        vec4 vec = loadu4(V); // { x, y, z, garbage }
        vec4 s0 = mul4(load4(val  ), vec);
        vec4 s1 = mul4(load4(val+4), vec);
        vec4 s2 = mul4(load4(val+8), vec);
        vec4 s3 = setzero4();
        s0 = add4(unpacklo4(s0, s1), unpackhi4(s0, s1));
        s2 = add4(unpacklo4(s2, s3), unpackhi4(s2, s3));
        return add4(catshift2(s0, s2), blend22(s0, s2));
    }

    /// multiplication by a 4-components vector: this * V
    const vec4 vecmul4_avx(vec4 const& vec) const
    {
        vec4 p = swap2f128(vec);
        vec4 l = blend22(vec, p);
        vec4 u = blend22(p, vec);
        vec4 x = mul4(load4(val   ), duplo4(l));
        vec4 y = mul4(load4(val+4 ), duphi4(l));
        x = fmadd4(load4(val+8 ), duplo4(u), x);
        y = fmadd4(load4(val+12), duphi4(u), y);
        return add4(x, y);
    }
    
    /// multiplication by a 4-components vector: transpose(this) * V
    const vec4 trans_vecmul4_avx(vec4 const& vec) const
    {
        vec4 s0 = mul4(load4(val   ), vec);
        vec4 s1 = mul4(load4(val+4 ), vec);
        vec4 s2 = mul4(load4(val+8 ), vec);
        vec4 s3 = mul4(load4(val+12), vec);
        s0 = add4(unpacklo4(s0, s1), unpackhi4(s0, s1));
        s2 = add4(unpacklo4(s2, s3), unpackhi4(s2, s3));
        return add4(catshift2(s0, s2), blend22(s0, s2));
    }

/*
    /// multiplication by a vector3: this * V
    const vec4 vecmul3_avx(real const* R) const
    {
        return vec4{val[0x0] * R[0] + val[0x4] * R[1] + val[0x8] * R[2],
                    val[0x1] * R[0] + val[0x5] * R[1] + val[0x9] * R[2],
                    val[0x2] * R[0] + val[0x6] * R[1] + val[0xA] * R[2], 0};
    }

    /// multiplication by a vector: transpose(M) * V
    const vec4 trans_vecmul3_avx(real const* R) const
    {
        return vec4{val[0x0] * R[0] + val[0x1] * R[1] + val[0x2] * R[2],
                    val[0x4] * R[0] + val[0x5] * R[1] + val[0x6] * R[2],
                    val[0x8] * R[0] + val[0x9] * R[1] + val[0xA] * R[2], 0};
    }

    /// multiplication by a vector: this * V
    const vec4 vecmul4_avx(real const* R) const
    {
        return vec4{val[0x0] * R[0] + val[0x4] * R[1] + val[0x8] * R[2] + val[0xC] * R[3],
                    val[0x1] * R[0] + val[0x5] * R[1] + val[0x9] * R[2] + val[0xD] * R[3],
                    val[0x2] * R[0] + val[0x6] * R[1] + val[0xA] * R[2] + val[0xE] * R[3],
                    val[0x3] * R[0] + val[0x7] * R[1] + val[0xB] * R[2] + val[0xF] * R[3]};
    }

    /// multiplication by a vector: transpose(M) * V
    const vec4 trans_vecmul4_avx(real const* R) const
    {
        return vec4{val[0x0] * R[0] + val[0x1] * R[1] + val[0x2] * R[2] + val[0x3] * R[3],
                    val[0x4] * R[0] + val[0x5] * R[1] + val[0x6] * R[2] + val[0x7] * R[3],
                    val[0x8] * R[0] + val[0x9] * R[1] + val[0xA] * R[2] + val[0xB] * R[3],
                    val[0xC] * R[0] + val[0xD] * R[1] + val[0xE] * R[2] + val[0xF] * R[3]};
    }
 */
#endif

    /// multiplication by a vector: this * V
    Vector3 vecmul3_(Vector3 const& V) const
    {
        return Vector3(val[0x0] * V.XX + val[0x4] * V.YY + val[0x8] * V.ZZ,
                       val[0x1] * V.XX + val[0x5] * V.YY + val[0x9] * V.ZZ,
                       val[0x2] * V.XX + val[0x6] * V.YY + val[0xA] * V.ZZ);
    }
    
    /// multiplication by a vector3: this * V
    Vector3 vecmul3_(real const* R) const
    {
        return Vector3(val[0x0] * R[0] + val[0x4] * R[1] + val[0x8] * R[2],
                       val[0x1] * R[0] + val[0x5] * R[1] + val[0x9] * R[2],
                       val[0x2] * R[0] + val[0x6] * R[1] + val[0xA] * R[2]);
    }

    /// multiplication by a vector: transpose(M) * V
    Vector3 trans_vecmul3_(Vector3 const& V) const
    {
        return Vector3(val[0x0] * V.XX + val[0x1] * V.YY + val[0x2] * V.ZZ,
                       val[0x4] * V.XX + val[0x5] * V.YY + val[0x6] * V.ZZ,
                       val[0x8] * V.XX + val[0x9] * V.YY + val[0xA] * V.ZZ);
    }
    
    /// multiplication by a vector: transpose(M) * V
    Vector3 trans_vecmul3_(real const* R) const
    {
        return Vector3(val[0x0] * R[0] + val[0x1] * R[1] + val[0x2] * R[2],
                       val[0x4] * R[0] + val[0x5] * R[1] + val[0x6] * R[2],
                       val[0x8] * R[0] + val[0x9] * R[1] + val[0xA] * R[2]);
    }

    /// multiplication by a vector: this * V
    Vector3 vecmul3(Vector3 const& vec) const
    {
#if MATRIX44_USES_AVX && VECTOR3_USES_AVX
        return vecmul3_avx(vec.vec);
#elif MATRIX44_USES_AVX
        return vecmul3_avx(vec.data());
#else
        return vecmul3_(vec);
#endif
    }

    /// multiplication by a vector: this * V
    Vector3 trans_vecmul3(Vector3 const& vec) const
    {
#if MATRIX44_USES_AVX
        return trans_vecmul3_avx(vec.data());
#else
        return trans_vecmul3_(vec);
#endif
    }
    
    /// vector multiplication
    friend Vector3 operator * (Matrix44 const& mat, Vector3 const& vec)
    {
        return mat.vecmul3(vec);
    }

    
    
    /// multiplication by a vector: this * V
    Vector4 vecmul4_(Vector4 const& V) const
    {
        return Vector4(val[0x0] * V.XX + val[0x4] * V.YY + val[0x8] * V.ZZ + val[0xC] * V.TT,
                       val[0x1] * V.XX + val[0x5] * V.YY + val[0x9] * V.ZZ + val[0xD] * V.TT,
                       val[0x2] * V.XX + val[0x6] * V.YY + val[0xA] * V.ZZ + val[0xE] * V.TT,
                       val[0x3] * V.XX + val[0x7] * V.YY + val[0xB] * V.ZZ + val[0xF] * V.TT);
    }
    
    /// multiplication by a vector: transpose(this) * V
    Vector4 vecmul4_(real const* R) const
    {
        return Vector4(val[0x0] * R[0] + val[0x4] * R[1] + val[0x8] * R[2] + val[0xC] * R[3],
                       val[0x1] * R[0] + val[0x5] * R[1] + val[0x9] * R[2] + val[0xD] * R[3],
                       val[0x2] * R[0] + val[0x6] * R[1] + val[0xA] * R[2] + val[0xE] * R[3],
                       val[0x3] * R[0] + val[0x7] * R[1] + val[0xB] * R[2] + val[0xF] * R[3]);
    }
    void vecmul4_add(real const* R, real* res) const
    {
        res[0x0]+=val[0x0] * R[0] + val[0x4] * R[1] + val[0x8] * R[2] + val[0xC] * R[3];
        res[0x1]+=val[0x1] * R[0] + val[0x5] * R[1] + val[0x9] * R[2] + val[0xD] * R[3];
        res[0x2]+=val[0x2] * R[0] + val[0x6] * R[1] + val[0xA] * R[2] + val[0xE] * R[3];
        res[0x3]+=val[0x3] * R[0] + val[0x7] * R[1] + val[0xB] * R[2] + val[0xF] * R[3];
    }
    
    /// multiplication by a vector: transpose(M) * V
    Vector4 trans_vecmul4_(Vector4 const& V) const
    {
        return Vector4(val[0x0] * V.XX + val[0x1] * V.YY + val[0x2] * V.ZZ + val[0x3] * V.TT,
                       val[0x4] * V.XX + val[0x5] * V.YY + val[0x6] * V.ZZ + val[0x7] * V.TT,
                       val[0x8] * V.XX + val[0x9] * V.YY + val[0xA] * V.ZZ + val[0xB] * V.TT,
                       val[0xC] * V.XX + val[0xD] * V.YY + val[0xE] * V.ZZ + val[0xF] * V.TT);
    }

    /// multiplication by a vector: transpose(M) * V
    Vector4 trans_vecmul4_(real const* R) const
    {
        return Vector4(val[0x0] * R[0] + val[0x1] * R[1] + val[0x2] * R[2] + val[0x3] * R[3],
                       val[0x4] * R[0] + val[0x5] * R[1] + val[0x6] * R[2] + val[0x7] * R[3],
                       val[0x8] * R[0] + val[0x9] * R[1] + val[0xA] * R[2] + val[0xB] * R[3],
                       val[0xC] * R[0] + val[0xD] * R[1] + val[0xE] * R[2] + val[0xF] * R[3]);
    }
    void trans_vecmul4_add(real const* R, real* res) const
    {
        res[0x0]+=val[0x0] * R[0] + val[0x1] * R[1] + val[0x2] * R[2] + val[0x3] * R[3];
        res[0x1]+=val[0x4] * R[0] + val[0x5] * R[1] + val[0x6] * R[2] + val[0x7] * R[3];
        res[0x2]+=val[0x8] * R[0] + val[0x9] * R[1] + val[0xA] * R[2] + val[0xB] * R[3];
        res[0x3]+=val[0xC] * R[0] + val[0xD] * R[1] + val[0xE] * R[2] + val[0xF] * R[3];
    }
    
   
             
     
            
     
    /*
         res2[0x0] += val[0x0]*R2[0];
         res1[0x0] += val[0x0]*R1[0];
         res1[0x1] += val[0x1]*R1[0];
         res2[0x0] += val[0x1]*R2[1];
         res2[0x0] += val[0x2]*R2[2];
         res1[0x2] += val[0x2]*R1[0];
         res2[0x0] += val[0x3]*R2[3];
         res1[0x3] += val[0x3]*R1[0];
         res1[0x0] += val[0x4]*R1[1];
         res2[0x1] += val[0x4]*R2[0];
         res2[0x1] += val[0x5]*R2[1];
         res1[0x1] += val[0x5]*R1[1];
         res1[0x2] += val[0x6]*R1[1];
         res2[0x1] += val[0x6]*R2[2];
         res1[0x3] += val[0x7]*R1[1];
         res2[0x1] += val[0x7]*R2[3];
            
         res1[0x0] += val[0x8]*  R1[2];
         res2[0x2] += val[0x8]*  R2[0];
         res1[0x1] += val[0x9] * R1[2];
         res2[0x2] += val[0x9] * R2[1];
         res2[0x2] += val[0xA] * R2[2];
         res1[0x2] += val[0xA]*  R1[2];
         res2[0x2] += val[0xB] * R2[3];
         res1[0x3] += val[0xB] * R1[2];
         res2[0x3] += val[0xC] * R2[0];
         res1[0x0] += val[0xC] * R1[3];
         res2[0x3] += val[0xD] * R2[1];
         res1[0x1] += val[0xD] * R1[3];
         res2[0x3] += val[0xE] * R2[2];
         res1[0x2] += val[0xE] * R1[3];
         res2[0x3] += val[0xF] * R2[3];
         res1[0x3] += val[0xF] * R1[3];
     */
    /*
    void bothvecmul(real const*R1, real const*R2, real*res, real*res2)
    {
        
        real coeff;
        real R1_val = R1[0];
        real temp_res2 = 0;
        real res1[4] = {0};
       
        
        coeff =val[0x0];
        temp_res2 += coeff*R2[0];
        res1[0x0] = coeff*R1_val;
        coeff =val[0x1];
        res1[0x1] += coeff*R1_val;
        temp_res2 += coeff*R2[1];
        coeff =val[0x2];
        temp_res2 += coeff*R2[2];
        res1[0x2] += coeff*R1_val;
        coeff =val[0x3];
        temp_res2 += coeff*R2[3];
        res1[0x3] += coeff*R1_val;
        coeff =val[0x4];
        
        res2[0x0]+= temp_res2;
        
        temp_res2 = 0;
        R1_val = R1[1];
        res1[0x0] += coeff*R1_val;
        temp_res2 += coeff*R2[0];
        
        coeff =val[0x5];
        temp_res2+= coeff*R2[1];
        res1[0x1] += coeff*R1_val;
        
        coeff =val[0x6];
        res1[0x2] += coeff*R1_val;
        temp_res2 += coeff*R2[2];
        
        coeff =val[0x7];
        res1[0x3] += coeff*R1_val;
        temp_res2 += coeff*R2[3];
        
        coeff =val[0x8];
        res2[0x1] +=temp_res2;
        temp_res2 = 0;
        R1_val = R1[2];
        res1[0x0] += coeff*  R1_val;
        temp_res2 += coeff*  R2[0];
        
        coeff =val[0x9];
        res1[0x1] += coeff * R1_val;
        temp_res2 += coeff * R2[1];
        
        coeff =val[0xA];
        temp_res2 += coeff * R2[2];
        res1[0x2] += coeff* R1_val;
        
        coeff =val[0xB];
        temp_res2 += coeff * R2[3];
        res1[0x3] += coeff * R1_val;
        
        coeff =val[0xC];
        R1_val = R1[3];
        res2[0x2]+= temp_res2;
        temp_res2 = 0;
        
        temp_res2 += coeff* R2[0];
        res1[0x0] += coeff* R1_val;
        coeff =val[0xD];
        
        temp_res2 += coeff * R2[1];
        res1[0x1] += coeff * R1_val;
        coeff =val[0xE];
        
        temp_res2 += coeff * R2[2];
        res1[0x2] += coeff * R1_val;
        
        coeff =val[0xF];
        temp_res2 += coeff* R2[3];
        
        res1[0x3] += coeff * R1_val;
        
        res2[0x3] += temp_res2;
        
        for(int i =0; i<4;i++)
        {
            res[i]+= res1[i];
        }
       
        
        
        
    }
     */
    void bothvecmul(real const* R1, real const*R2, real*res, real*res3)
    {
        real y10 = 0;
        real y11 = 0;
        real y12=0;
        real y13=0;
        real y20 = 0;
        real y21 = 0;
        real y22 = 0;
        real y23=0;
            real r10 = R1[0];
                real r11 = R1[1];
                real r12 = R1[2];
                real r13 = R1[3];
                real r20 = R2[0];
                real r21 = R2[1];
                real r22 = R2[2];
                real r23 = R2[3];
                {real & c = val[0x0];

                y10 += c*r10;
                y20 += c*r20;
                }
                {
                real & c = val[0x1];
                y11 += c * r10;
                y20 += c * r21;
                }
                {
                real & c = val[0x4];
                y10 += c * r11;
                y21 += c * r20;
                }
                {
                real & c = val[0x2];
                y12 += c * r10;
                y20 += c * r22;
                }
                {
                real& c = val[0x8];
                y10 += c * r12;
                y22 += c* r20;
                }
                {
                real & c = val[0x3];
                y13 += c * r10;
                y20 += c * r23;
                }
                {
                real & c  = val[0xC];
                y10 += c *r13;
                y23 += c * r20;
                }
                {
                real & c = val[0x5];
                y11+= c*r11;
                y21 += c *r21;
                }
                {
                real & c=  val[0x6];
                y12 += c*r11;
                y21 += c*r22;
                }
                {
                real & c = val[0x9];
                y11 += c *r12;
                y22 += c *r21;
                }
                {
                real &c  = val[0x7];
                y13 += c *r11;
                y21 += c * r23;
                }
                {
                real & c  = val[0xD];
                y11 +=c*r13;
                y23 += c *r21;
                }
                {
                real & c  = val[0xA];
                y12 += c * r12;
                y22 += c* r22;
                }
                {
                real & c = val[0xE];
                y12 += c * r13;
                y23 += c * r22;
                }
                {
                real & c = val[0xB];
                y13 += c * r12;
                y22 += c *r23;
                }
                {
                real & c = val[0xF];
                y13 += c * r13;
                y23 += c *r23;
                }
        
                res[0] += y10;
                res[1] += y11;
                res[2] += y12;
                res[3] += y13;
                res3[0] += y20;
                res3[1] += y21;
                res3[2] += y22;
                res3[3] += y23;
                
            
    }
    /// multiplication by a vector: this * V
    Vector4 vecmul(Vector4 const& vec) const
    {
#if MATRIX44_USES_AVX && VECTOR3_USES_AVX
        return vecmul4_avx(vec.vec);
#elif MATRIX44_USES_AVX
        return vecmul4_avx(vec.data());
#else
        return vecmul4_(vec);
#endif
    }

    /// multiplication by a vector: this * V
    Vector4 trans_vecmul(Vector4 const& vec) const
    {
#if MATRIX44_USES_AVX && VECTOR3_USES_AVX
        return trans_vecmul4_avx(vec.vec);
#elif MATRIX44_USES_AVX
        return trans_vecmul4_avx(vec.data());
#else
        return trans_vecmul4_(vec);
#endif
    }

    /// vector multiplication
    friend Vector4 operator * (Matrix44 const& mat, Vector4 const& vec)
    {
        return mat.vecmul(vec);
    }
    

    /// multiplication by another matrix: @returns this * M
    const Matrix44 mul(Matrix44 const& M) const
    {
        Matrix44 res;
        res[0x0] = val[0x0] * M[0x0] + val[0x4] * M[0x1] + val[0x8] * M[0x2] + val[0xC] * M[0x3];
        res[0x1] = val[0x1] * M[0x0] + val[0x5] * M[0x1] + val[0x9] * M[0x2] + val[0xD] * M[0x3];
        res[0x2] = val[0x2] * M[0x0] + val[0x6] * M[0x1] + val[0xA] * M[0x2] + val[0xE] * M[0x3];
        res[0x3] = val[0x3] * M[0x0] + val[0x7] * M[0x1] + val[0xB] * M[0x2] + val[0xF] * M[0x3];

        res[0x4] = val[0x0] * M[0x4] + val[0x4] * M[0x5] + val[0x8] * M[0x6] + val[0xC] * M[0x7];
        res[0x5] = val[0x1] * M[0x4] + val[0x5] * M[0x5] + val[0x9] * M[0x6] + val[0xD] * M[0x7];
        res[0x6] = val[0x2] * M[0x4] + val[0x6] * M[0x5] + val[0xA] * M[0x6] + val[0xE] * M[0x7];
        res[0x7] = val[0x3] * M[0x4] + val[0x7] * M[0x5] + val[0xB] * M[0x6] + val[0xF] * M[0x7];

        res[0x8] = val[0x0] * M[0x8] + val[0x4] * M[0x9] + val[0x8] * M[0xA] + val[0xC] * M[0xB];
        res[0x9] = val[0x1] * M[0x8] + val[0x5] * M[0x9] + val[0x9] * M[0xA] + val[0xD] * M[0xB];
        res[0xA] = val[0x2] * M[0x8] + val[0x6] * M[0x9] + val[0xA] * M[0xA] + val[0xE] * M[0xB];
        res[0xB] = val[0x3] * M[0x8] + val[0x7] * M[0x9] + val[0xB] * M[0xA] + val[0xF] * M[0xB];

        res[0xC] = val[0x0] * M[0xC] + val[0x4] * M[0xD] + val[0x8] * M[0xE] + val[0xC] * M[0xF];
        res[0xD] = val[0x1] * M[0xC] + val[0x5] * M[0xD] + val[0x9] * M[0xE] + val[0xD] * M[0xF];
        res[0xE] = val[0x2] * M[0xC] + val[0x6] * M[0xD] + val[0xA] * M[0xE] + val[0xE] * M[0xF];
        res[0xF] = val[0x3] * M[0xC] + val[0x7] * M[0xD] + val[0xB] * M[0xE] + val[0xF] * M[0xF];
        return res;
    }
    
    /// multiplication by matrix
    friend Matrix44 operator * (Matrix44 const& mat, Matrix44 const& mut)
    {
        return mat.mul(mut);
    }

    /// multiplication by another matrix: @returns transpose(this) * M
    const Matrix44 trans_mul(Matrix44 const& M) const
    {
        ABORT_NOW("unfinished");
    }
    
    /// add full matrix: this <- this + M
    void add_full(Matrix44 const& M)
    {
        real const* src = M.val;
        for ( index u = 0; u < 16; ++u )
            val[u] += src[u];
    }
    
    /// add full matrix: this <- this + alpha * M
    void add_full(const real alpha, Matrix44 const& M)
    {
        real const* src = M.val;
        for ( index u = 0; u < 16; ++u )
            val[u] += alpha * src[u];
    }
    
    /// sub full matrix: this <- this - M
    void sub_full(Matrix44 const& M)
    {
        real const* src = M.val;
        for ( index u = 0; u < 16; ++u )
            val[u] -= src[u];
    }
    
    /// subtract full matrix: this <- this - alpha * M
    void sub_full(const real alpha, Matrix44 const& M)
    {
        real const* src = M.val;
        for ( index u = 0; u < 16; ++u )
            val[u] -= alpha * src[u];
    }

    /// add lower triangle of matrix including diagonal: this <- this + M
    void add_half(Matrix44 const& M)
    {
        real const* src = M.val;
#if ( 1 )
        for ( index u = 0; u < 16; ++u )
            val[u] += src[u];
#else
        for ( index x = 0; x < 4; ++x )
        for ( index y = x; y < 4; ++y )
            val[y+4*x] += src[y+4*x];
#endif
    }
    
    /// add lower triangle of matrix including diagonal: this <- this + alpha * M
    void add_half(const real alpha, Matrix44 const& M)
    {
        real const* src = M.val;
#if ( 1 )
        for ( index u = 0; u < 16; ++u )
            val[u] += alpha * src[u];
#else
        for ( index x = 0; x < 4; ++x )
        for ( index y = x; y < 4; ++y )
            val[y+4*x] += alpha * src[y+4*x];
#endif
    }
    
    /// add alpha to diagonal
    void add_diag(real alpha)
    {
        val[0x0] += alpha;
        val[0x5] += alpha;
        val[0xA] += alpha;
        val[0xF] += alpha;
    }
    
    /// add -alpha to diagonal
    void sub_diag(real alpha)
    {
        val[0x0] -= alpha;
        val[0x5] -= alpha;
        val[0xA] -= alpha;
        val[0xF] -= alpha;
    }

    /// subtract lower triangle of matrix including diagonal: this <- this - M
    void sub_half(Matrix44 const& M)
    {
        real const* src = M.val;
#if ( 1 )
        for ( index u = 0; u < 16; ++u )
            val[u] -= src[u];
#else
        for ( index x = 0; x < 4; ++x )
        for ( index y = x; y < 4; ++y )
            val[y+4*x] -= src[y+4*x];
#endif
    }

    
    /// add all elements of block 'S' to array 'M'
    void addto(real * M, size_t ldd) const
    {
        for ( index x = 0; x < 4; ++x )
        for ( index y = 0; y < 4; ++y )
            M[y+ldd*x] = val[y+4*x];
    }
    
    /// add scaled elements of block 'S' to array 'M'
    void addto(real * M, size_t ldd, real alpha) const
    {
        for ( index x = 0; x < 4; ++x )
        for ( index y = 0; y < 4; ++y )
            M[y+ldd*x] = alpha * val[y+4*x];
    }

    /// add lower elements of this block to lower triangle of 'M'
    void addto_lower(real * M, size_t ldd) const
    {
        for ( index x = 0; x < 4; ++x )
        for ( index y = x; y < 4; ++y )
            M[y+ldd*x] = val[y+4*x];
    }
    
    /// add scaled lower elements of this block to lower triangle of 'M'
    void addto_lower(real * M, size_t ldd, real alpha) const
    {
        for ( index x = 0; x < 4; ++x )
        for ( index y = x; y < 4; ++y )
            M[y+ldd*x] = alpha * val[y+4*x];
    }

    /// add all elements of this block to 'M', with transposition
    void addto_trans(real * M, size_t ldd) const
    {
        for ( index x = 0; x < 4; ++x )
        for ( index y = 0; y < 4; ++y )
            M[x+ldd*y] = val[y+4*x];
    }
    
    /// add lower elements of this block to both upper and lower triangles of 'M'
    void addto_symm(real * M, size_t ldd) const
    {
        for ( index x = 0; x < 4; ++x )
        {
            M[x+ldd*x] = val[x+4*x];
            for ( index y = x+1; y < 4; ++y )
            {
                M[y+ldd*x] = val[y+4*x];
                M[x+ldd*y] = val[y+4*x];
            }
        }
    }

#pragma mark -
    
    /// return diagonal Matrix from diagonal terms
    static Matrix44 diagonal(real a, real b, real c, real d)
    {
        return Matrix44(a, 0, 0, 0, 0, b, 0, 0, 0, 0, c, 0, 0, 0, 0, d);
    }
    
    /// identity matrix
    static Matrix44 one()
    {
        return Matrix44(0, 1);
    }

    /// construct Matrix from coordinates (column-major)
    static Matrix44 symmetric(real a, real b, real c, real d,
                              real e, real f, real g, real h,
                              real i, real j )
    {
        return Matrix44(a, b, c, d, b, e, f, g, c, f, h, i, d, g, i, j);
    }

    /// return a symmetric matrix: [ dir (x) transpose(dir) ]
    static Matrix44 outerProduct(Vector4 const& V)
    {
        return symmetric(V[0]*V[0], V[1]*V[0], V[2]*V[0], V[3]*V[0],
                         V[1]*V[1], V[2]*V[1], V[3]*V[1],
                         V[2]*V[2], V[3]*V[2],
                         V[3]*V[3] );
    }
    
    /// return a symmetric matrix: alpha * [ dir (x) transpose(dir) ]
    static Matrix44 outerProduct(Vector4 const& V, real alpha)
    {
        real X = V[0] * alpha;
        real Y = V[1] * alpha;
        real Z = V[2] * alpha;
        real T = V[3] * alpha;
        return symmetric(V[0]*X, V[1]*X, V[2]*X, V[3]*X,
                         V[1]*Y, V[2]*Y, V[3]*Y,
                         V[2]*Z, V[3]*Z,
                         V[3]*T );
    }
    
    /// return outer product: [ dir (x) transpose(vec) ]
    static Matrix44 outerProduct(const real D[], const real V[])
    {
#if MATRIX44_USES_AVX
        Matrix44 res;
        vec4 s = load4(V);
        vec4 p = swap2f128(s);
        vec4 l = blend22(s, p);
        vec4 u = blend22(p, s);
        vec4 d = load4(D);
        store4(res.val   , mul4(d, duplo4(l)));
        store4(res.val+4 , mul4(d, duphi4(l)));
        store4(res.val+8 , mul4(d, duplo4(u)));
        store4(res.val+12, mul4(d, duphi4(u)));
        return res;
#else
        return Matrix44(D[0]*V[0], D[1]*V[0], D[2]*V[0], D[3]*V[0],
                        D[0]*V[1], D[1]*V[1], D[2]*V[1], D[3]*V[1],
                        D[0]*V[2], D[1]*V[2], D[2]*V[2], D[3]*V[2],
                        D[0]*V[3], D[1]*V[3], D[2]*V[3], D[3]*V[3] );
#endif
    }
    
    /// return [ dir (x) transpose(vec) + vec (x) transpose(dir) ]
    static Matrix44 symmetricOuterProduct(const real D[], const real V[])
    {
        real xx = D[0] * V[0];
        real yy = D[1] * V[1];
        real zz = D[2] * V[2];
        real tt = D[3] * V[3];
        return symmetric(xx+xx, D[1]*V[0] + D[0]*V[1], D[2]*V[0] + D[0]*V[2], D[3]*V[0] + D[0]*V[3],
                         yy+yy, D[2]*V[1] + D[1]*V[2], D[3]*V[1] + D[1]*V[3],
                         zz+zz, D[3]*V[2] + D[2]*V[3],
                         tt+tt);
    }
 
    /// return symmetric matrix block :  dia * I + [ dir (x) dir ] * len
    static Matrix44 offsetOuterProduct(const real dia, Vector4 const& dir, const real len)
    {
        real xl = dir.XX * len;
        real yl = dir.YY * len;
        real zl = dir.ZZ * len;
        real tl = dir.TT * len;
        return symmetric(xl * dir.XX + dia, yl * dir.XX, zl * dir.XX, tl * dir.XX,
                         yl * dir.YY + dia, zl * dir.YY, tl * dir.YY,
                         zl * dir.ZZ + dia, tl * dir.ZZ,
                         tl * dir.TT + dia);
    }
    
    /// return symmetric matrix block :  dia * I + [ dir (x) dir ] * len
    static Matrix44 offsetOuterProduct(const real dia, const real D[], const real len)
    {
        real xl = D[0] * len;
        real yl = D[1] * len;
        real zl = D[2] * len;
        real tl = D[3] * len;
        return symmetric(xl * D[0] + dia, yl * D[0], zl * D[0], tl * D[0],
                         yl * D[1] + dia, zl * D[1], tl * D[1],
                         zl * D[2] + dia, tl * D[2],
                         tl * D[3] + dia);
    }
};


/// output operator
inline std::ostream& operator << (std::ostream& os, Matrix44 const& arg)
{
    std::ios::fmtflags fgs = os.flags();
    arg.print(os);
    os.setf(fgs);
    return os;
}

#endif

