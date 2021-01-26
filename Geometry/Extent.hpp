/////////////////////////////////////////////////////////////////////
// Objects and methods for representing and managing object        //
// extents.                                                        //
/////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////
// -- HISTORY ---------------------------------------------------- //
// 01/22/2021 - Brennan Young                                      //
// - created                                                       //
// 01/25/2021 - Brennan Young                                      //
// - default constructor now initializes up to dimension n.        //
/////////////////////////////////////////////////////////////////////

#ifndef YOUNG_GEOMETRY_EXTENT_20210122
#define YOUNG_GEOMETRY_EXTENT_20210122

#include <cstddef> // size_t


class Extent {
private:
    size_t n;     // number of dimensions
    double *a, *b; // minimum, maximum bounds
    
public:
    // constructor
    Extent ( size_t N=1 ) : n(N)
    {
        if ( n < 1 ) n = 1;
        a = new double [n];
        b = new double [n];
        for ( size_t i = 0; i < n; ++i ) a[i] = b[i] = 0.0;
    }
    
    // 1D constructor
    Extent ( double x0, double x1, double buff=0.0 ) : n(1)
    {
        a = new double [n];
        b = new double [n];
        if ( x0 < x1 ) { a[0] = x0;  b[0] = x1; }
        else           { a[0] = x1;  b[0] = x0; }
        for ( size_t i = 0; i < n; ++i ) {
            a[i] -= buff;
            b[i] += buff;
        }
    }
    
    // 2D constructor
    Extent ( double x0, double y0, double x1, double y1,
        double buff=0.0 )
    : n(2)
    {
        a = new double [n];
        b = new double [n];
        if ( x0 < x1 ) { a[0] = x0;  b[0] = x1; }
        else           { a[0] = x1;  b[0] = x0; }
        if ( y0 < y1 ) { a[1] = y0;  b[1] = y1; }
        else           { a[1] = y1;  b[1] = y0; }
        for ( size_t i = 0; i < n; ++i ) {
            a[i] -= buff;
            b[i] += buff;
        }
    }
    
    // 3D constructor
    Extent ( double x0, double y0, double z0,
        double x1, double y1, double z1, double buff=0.0 )
    : n(3)
    {
        a = new double [n];
        b = new double [n];
        if ( x0 < x1 ) { a[0] = x0;  b[0] = x1; }
        else           { a[0] = x1;  b[0] = x0; }
        if ( y0 < y1 ) { a[1] = y0;  b[1] = y1; }
        else           { a[1] = y1;  b[1] = y0; }
        if ( z0 < z1 ) { a[2] = z0;  b[2] = z1; }
        else           { a[2] = z1;  b[2] = z0; }
        for ( size_t i = 0; i < n; ++i ) {
            a[i] -= buff;
            b[i] += buff;
        }
    }
    
    // copy constructor
    Extent ( const Extent& E ) : n(E.n)
    {
        a = new double [n];
        b = new double [n];
        for ( size_t i = 0; i < n; ++i ) {
            a[i] = E.a[i];
            b[i] = E.b[i];
        }
    }
    
    // destructor
    ~Extent ()
    {
        delete[] a;
        delete[] b;
    }
    
    // operators
    Extent& operator= ( const Extent& E )
    {
        if ( &E == this ) return *this;
        
        delete[] a;
        delete[] b;
        n = E.n;
        a = new double [n];
        b = new double [n];
        for ( size_t i = 0; i < n; ++i ) {
            a[i] = E.a[i];
            b[i] = E.b[i];
        }
        
        return *this;
    }
    
    // less-than comparison
    bool operator< ( const Extent& E ) const
    {
        // compare size
        if ( n < E.n ) return true;
        if ( n > E.n ) return false;
        
        // compare minimum bound
        for ( size_t i = 0; i < n; ++i ) {
            if ( a[i] < E.a[i] ) return true;
            if ( a[i] > E.a[i] ) return false;
        }
        
        // compare maximum bound
        for ( size_t i = 0; i < n; ++i ) {
            if ( b[i] < E.b[i] ) return true;
            if ( b[i] > E.b[i] ) return false;
        }
        
        return false; // same
    }
    
    // greater-than comparison
    bool operator> ( const Extent& E ) const
    {
        // compare size
        if ( n < E.n ) return false;
        if ( n > E.n ) return true;
        
        // compare minimum bound
        for ( size_t i = 0; i < n; ++i ) {
            if ( a[i] < E.a[i] ) return false;
            if ( a[i] > E.a[i] ) return true;
        }
        
        // compare maximum bound
        for ( size_t i = 0; i < n; ++i ) {
            if ( b[i] < E.b[i] ) return false;
            if ( b[i] > E.b[i] ) return true;
        }
        
        return false; // same
    }
    
    // equal-to comparison
    bool operator== ( const Extent& E ) const
    {
        size_t i = 0;
        for ( ; i < n && i < E.n; ++i )
            if ( a[i] != E.a[i] || b[i] != E.b[i] ) return false;
        return !(i < n || i < E.n);
    }
    bool operator!= ( const Extent& E ) const
    {
        size_t i = 0;
        for ( ; i < n && i < E.n; ++i )
            if ( a[i] != E.a[i] || b[i] != E.b[i] ) return true;
        return i < n || i < E.n;
    }
    
    // get number of dimensions
    size_t size () const { return n; }
    
    // element access
    double& min ( size_t i )       { return a[i]; }
    double  min ( size_t i ) const { return a[i]; }
    double& max ( size_t i )       { return b[i]; }
    double  max ( size_t i ) const { return b[i]; }
    
    double& xmin () const { return a[0]; }
    double  xmin ()       { return a[0]; }
    double& xmax () const { return b[0]; }
    double  xmax ()       { return b[0]; }
    double& ymin () const { return a[1]; }
    double  ymin ()       { return a[1]; }
    double& ymax () const { return b[1]; }
    double  ymax ()       { return b[1]; }
    double& zmin () const { return a[2]; }
    double  zmin ()       { return a[2]; }
    double& zmax () const { return b[2]; }
    double  zmax ()       { return b[2]; }
    
    // get the length in a dimension
    double len ( size_t i ) const { return b[i] - a[i]; }
    double xlen () const { return b[0] - a[0]; }
    double ylen () const { return b[1] - a[1]; }
    double zlen () const { return b[2] - a[2]; }
    
    // get the (hyper-volume) of the extent
    double volume () const
    {
        double sum = 1.0;
        for ( size_t i = 0; i < n; ++i ) sum *= len(i);
        return sum;
    }
    
    // expand
    Extent& add ( const Extent& E )
    {
        for ( size_t i = 0; i < n && i < E.n; ++i ) {
            if ( E.a[i] < a[i] ) a[i] = E.a[i];
            if ( E.b[i] > b[i] ) b[i] = E.b[i];
        }
        return *this;
    }
    Extent add ( const Extent& E ) const
    {
        Extent ext = *this;
        return ext.add(E);
    }
    Extent operator+ ( const Extent& E ) const { return add(E); }
    Extent& operator+= ( const Extent& E ) { return add(E); }
    
    // reduce this extent to the intersection with another extent
    Extent& clip ( const Extent& E )
    {
        for ( size_t i = 0; i < n && i < E.n; ++i ) {
            if      ( b[i] < E.a[i] ) a[i] = b[i] = 0.0;
            else if ( E.b[i] < a[i] ) a[i] = b[i] = 0.0;
            else {
                a[i] = a[i] < E.a[i] ? E.a[i] : a[i];
                b[i] = b[i] > E.b[i] ? E.b[i] : b[i];
            }
        }
        return *this;
    }
    Extent& operator-= ( const Extent& E ) { return clip(E); }
    
    // get the bounds of the intersection of two extents
    Extent intersect ( const Extent& E ) const
    {
        Extent ext = *this;
        return ext.clip(E);
    }
    Extent operator- ( const Extent& E ) const
    {
        return intersect(E);
    }
    
    // determine if overlaps with another extent's bounds
    bool overlap ( const Extent& E ) const
    {
        return intersect(E).volume() > 0.0;
    }
}; // Extent


#endif // YOUNG_GEOMETRY_EXTENT_20210122