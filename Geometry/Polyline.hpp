/////////////////////////////////////////////////////////////////////
// Objects and methods for managing geometric polylines (2d).      //
/////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////
// -- HISTORY ---------------------------------------------------- //
// 12/29/2020 - Brennan Young                                      //
// - created                                                       //
// 01/13/2020 - Brennan Young                                      //
// - added constructor with Polygon2d argument.                    //
// - added extractChain method.                                    //
// - added removeChain method.                                     //
/////////////////////////////////////////////////////////////////////

#ifndef YOUNG_GEOMETRY_POLYLINE_20210111
#define YOUNG_GEOMETRY_POLYLINE_20210111

#include <vector>      // std::vector, size_t
#include "Polygon.hpp" // Point2d, Line2d, Polygon2d


class Polyline2d {
public:
    std::vector<Point2d> v; // vertices
    std::vector<size_t> s;  // indexes in v that start a new chain
    Point2d bMin, bMax;     // bounds
    
    // constructors, destructor
    Polyline2d (
        const std::vector<Point2d>& V=std::vector<Point2d>() ) : v(V)
    {
        computeBounds();
    }
    Polyline2d ( const Polygon2d& G ) : v(G.v), bMin(G.bMin),
        bMax(G.bMax)
    {
        size_t i0, i1=0;
        
        // get the beginning of each chain
        for ( size_t i = 0; i < G.size(); ++i ) {
            // start of polygon
            if ( i == i1 ) {
                i0 = i;
                s.push_back(i);
            }
            
            // end of polygon
            else if ( G[i] == G[i0] ) {
                i1 = i+1;
                continue;
            }
        }
    }
    Polyline2d ( const Polyline2d& L ) : v(L.v), s(L.s),
        bMin(L.bMin), bMax(L.bMax) {}
    ~Polyline2d () {}
    
    // operators
    Polyline2d& operator= ( const Polyline2d& L )
    {
        if ( &L == this ) return *this;
        v = L.v;
        s = L.s;
        bMin = L.bMin;
        bMax = L.bMax;
        return *this;
    }
    Point2d& operator[] ( size_t i )
    {
        return v[i];
    }
    const Point2d& operator[] ( size_t i ) const
    {
        return v[i];
    }
    bool operator< ( const Polyline2d& L ) const
    {
        // sort by coordinates first
        if ( bMin < L.bMin ) return true;
        if ( bMin > L.bMin ) return false;
        if ( bMax < L.bMax ) return true;
        if ( bMax > L.bMax ) return false;
        // then sort by number of chains
        if ( s.size() < L.s.size() ) return true;
        if ( s.size() > L.s.size() ) return false;
        // then sort by number of vertices
        if ( v.size() < L.v.size() ) return true;
        if ( v.size() > L.v.size() ) return false;
        // then sort by the vertices themselves
        for ( size_t i = 0; i < v.size(); ++i ) {
            if ( v[i] < L.v[i] ) return true;
            if ( v[i] > L.v[i] ) return false;
        }
        return false;
    }
    bool operator== ( const Polyline2d& L ) const
    {
        size_t i = 0;
        
        // compare chains
        for ( ; i < s.size() && i < L.s.size(); ++i )
            if ( s[i] != L.s[i] ) return false;
        if ( i < s.size() || i < L.s.size() ) return false;
        
        // compare vertices
        for ( i = 0; i < v.size() && i < L.v.size(); ++i )
            if ( v[i] != L.v[i] ) return false;
        return !(i < v.size() || i < L.v.size());
    }
    bool operator!= ( const Polyline2d& L ) const
    {
        size_t i = 0;
        
        // compare chains
        for ( ; i < s.size() && i < L.s.size(); ++i )
            if ( s[i] != L.s[i] ) return true;
        if ( i < s.size() || i < L.s.size() ) return true;
        
        // compare vertices
        for ( i = 0; i < v.size() && i < L.v.size(); ++i )
            if ( v[i] != L.v[i] ) return true;
        return i < v.size() || i < L.v.size();
    }
    
    // getters
    size_t size () const
    {
        return v.size();
    }
    
    size_t chains () const
    {
        return s.size();
    }
    
    // basic operations
    
    // create a polyline of chain i
    Polyline2d extractChain ( size_t i ) const
    {
        if ( i >= s.size() ) return Polyline2d();
        
        size_t jend = (i+1 == s.size() ? v.size() : s[i+1]);
        
        Polyline2d out;
        out.v.reserve(jend - s[i]);
        
        for ( size_t j = s[i]; j < jend; ++j )
            out.v.push_back(v[j]);
        out.computeBounds();
        
        return out;
    }
    
    // delete chain i from polyline
    void removeChain ( size_t i )
    {
        if ( i >= s.size() ) return;
        
        size_t jend = (i+1 == s.size() ? v.size() : s[i+1]);
        
        std::vector<Point2d> vnew;
        std::vector<Point2d> snew;
        vnew.reserve(v.size() - (jend - s[i]));
        snew.reserve(s.size() - 1);
        
        // create new vector of chains to keep
        for ( size_t I = 0; I < s.size(); ++I ) {
            if ( I == i ) continue;
            size_t jend = (I+1 == s.size() ? v.size() : s[I+1]);
            
            // add chain
            snew.push_back(I);
            for ( size_t j = s[I]; j < jend; ++j )
                vnew.push_back(v[j]);
        }
        
        // replace old vectors
        v = vnew;
        s = snew;
        computeBounds();
    }
    
    // compute minimum and maximum bounding points; ensure chains
    // are represented in data structure
    void computeBounds ()
    {
        if ( v.size() > 0 ) {
            if ( s.size() == 0 ) s.push_back(0);
            bMin = bMax = v[0];
            for ( size_t i = 1; i < v.size(); ++i ) {
                if ( v[i].x < bMin.x ) bMin.x = v[i].x;
                else if ( v[i].x > bMax.x ) bMax.x = v[i].x;
                if ( v[i].y < bMin.y ) bMin.y = v[i].y;
                else if ( v[i].y > bMax.y ) bMax.y = v[i].y;
            }
        }
        else {
            s.clear();
            bMin = bMax = Point2d();
        }
    }
}; // Polyline2d


#endif // YOUNG_GEOMETRY_POLYLINE_20210111