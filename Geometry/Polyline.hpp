/////////////////////////////////////////////////////////////////////
// Objects and methods for managing geometric polylines (2d).      //
//                                                                 //
// Polylines may represent multiple chains of vertices -- the      //
// beginning of each chain *must* be represented in the 's' member //
// in order (s[i] < s[i+1]).                                       //
//                                                                 //
// Polyline chains may loop on themselves. These *must* be         //
// represented in the 'p' member, which must be of the same length //
// and in the same order as 's'. Loops *must* have equivalent      //
// starting and ending vertices.                                   //
/////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////
// -- HISTORY ---------------------------------------------------- //
// 12/29/2020 - Brennan Young                                      //
// - created                                                       //
// 01/13/2020 - Brennan Young                                      //
// - added constructor with Polygon2d argument.                    //
// - added extractChain method.                                    //
// - added removeChain method.                                     //
// 01/20/2021 - Brennan Young                                      //
// - clean up minor errors.                                        //
// - added operators + and += for appending polylines.             //
// - added loop, to indicate whether or not a chain loops (i.e.,   //
//   has no "first" or "last" vertex).                             //
// 01/21/2021 - Brennan Young                                      //
// - overloaded addChain for adding a new chain at the tail of the //
//   polyline.                                                     //
// - added append for appending vertices to an existing chain.     //
// 01/21/2021 - Brennan Young                                      //
// - better attention has been paid to updating bounds.            //
// 01/22/2021 - Brennan Young                                      //
// - replaced Point2d bMin and bMax objects with Extent ext.       //
// 01/25/2021 - Brennan Young                                      //
// - made v, s, loop, and ext private, in order to protect vertex  //
//   and chain management from misalignment. Added relevant        //
//   and setters and getters.                                      //
// - added private Chain struct to help manage chains.             //
// - renamed chains() to nchains().                                //
// - renamed addChain(size_t,bool) to addChainBreak.               //
// - added removeChainBreak(size_t).                               //
// - other minor adjustments.                                      //
/////////////////////////////////////////////////////////////////////

#ifndef YOUNG_GEOMETRY_POLYLINE_20210111
#define YOUNG_GEOMETRY_POLYLINE_20210111

#include <vector>      // std::vector, size_t
#include "Polygon.hpp" // Extent, Point2d, Line2d, Polygon2d


class Polyline2d {
private:
    // structure
    struct Chain {
        size_t start;   // polyline index of first vertex in chain
        bool   loop;    // true if the chain loops
        Chain (size_t s=0, bool p=false) : start(s), loop(p) {}
    }; // Polyline2d::Chain
    
    std::vector<Point2d> v; // vertices
    std::vector<Chain> s;   // chain breaks
    Extent ext;             // bounds

public:
    // constructor
    Polyline2d (
        const std::vector<Point2d>& V=std::vector<Point2d>(),
        bool c=false )
    : v(V), ext(Extent(2))
    {
        s.push_back(Chain(0, c));
        computeBounds();
    }
    
    // construct from polygon
    Polyline2d ( const Polygon2d& G )
    : v(G.v), ext(G.ext)
    {
        size_t i0, i1=0;
        
        // get the beginning of each chain
        for ( size_t i = 0; i < G.size(); ++i ) {
            // start of polygon
            if ( i == i1 ) {
                i0 = i;
                s.push_back(Chain(i, true));
            }
            
            // end of polygon
            else if ( G[i] == G[i0] ) {
                i1 = i+1;
                continue;
            }
        }
    }
    
    // copy constructor
    Polyline2d ( const Polyline2d& L ) : v(L.v), s(L.s), ext(L.ext)
    {}
    
    // destructor
    ~Polyline2d () {}
    
    // operators
    
    // assignment
    Polyline2d& operator= ( const Polyline2d& L )
    {
        if ( &L == this ) return *this;
        v = L.v;
        s = L.s;
        ext = L.ext;
        return *this;
    }
    
    // element access
    Point2d& operator[] ( size_t i ) { return v[i]; }
    const Point2d& operator[] ( size_t i ) const { return v[i]; }
    
    // less-than comparison
    bool operator< ( const Polyline2d& L ) const
    {
        // sort by coordinates first
        if ( ext < L.ext ) return true;
        if ( ext > L.ext ) return false;
        
        // then sort by number of chains
        if ( s.size() < L.s.size() ) return true;
        if ( s.size() > L.s.size() ) return false;
        
        // then sort by number of vertices
        if ( v.size() < L.v.size() ) return true;
        if ( v.size() > L.v.size() ) return false;
        
        // then sort by chain places
        for ( size_t i = 0; i < s.size(); ++i ) {
            if ( s[i].start < L.s[i].start ) return true;
            if ( s[i].start > L.s[i].start ) return false;
        }
        
        // then sort by chain looping
        for ( size_t i = 0; i < s.size(); ++i ) {
            if ( s[i].loop < L.s[i].loop ) return true;
            if ( s[i].loop > L.s[i].loop ) return false;
        }
        
        // then sort by the vertices themselves
        for ( size_t i = 0; i < v.size(); ++i ) {
            if ( v[i] < L.v[i] ) return true;
            if ( v[i] > L.v[i] ) return false;
        }
        
        return false; // same
    }
    
    // equal to comparison
    bool operator== ( const Polyline2d& L ) const
    {
        size_t i = 0;
        
        // compare chains
        for ( ; i < s.size() && i < L.s.size(); ++i ) {
            if ( s[i].start != L.s[i].start
                    || s[i].loop != L.s[i].loop )
                return false;
        }
        if ( i < s.size() || i < L.s.size() ) return false;
        
        // compare vertices
        for ( i = 0; i < v.size() && i < L.v.size(); ++i )
            if ( v[i] != L.v[i] ) return false;
        return !(i < v.size() || i < L.v.size());
    }
    
    // not equal to comparison
    bool operator!= ( const Polyline2d& L ) const
    {
        size_t i = 0;
        
        // compare chains
        for ( ; i < s.size() && i < L.s.size(); ++i )
            if ( s[i].start != L.s[i].start ) return true;
        if ( i < s.size() || i < L.s.size() ) return true;
        
        // compare vertices
        for ( i = 0; i < v.size() && i < L.v.size(); ++i )
            if ( v[i] != L.v[i] ) return true;
        return i < v.size() || i < L.v.size();
    }
    
    // append polyline
    Polyline2d operator+ ( const Polyline2d& L ) const
    {
        Polyline2d P (*this);
        
        // add chains
        size_t v0 = v.size();
        P.v.insert(P.v.end(), L.v.begin(), L.v.end());
        P.s.reserve(P.s.size() + L.s.size());
        for ( size_t i = 0; i < L.s.size(); ++i ) {
            P.s.push_back(Chain(v0 + L.s[i].start, L.s[i].loop));
        }
        
        // update bounds
        ext.add(L.ext);
        
        return P;
    }
    Polyline2d& operator+= ( const Polyline2d& L )
    {
        // add chains
        size_t v0 = v.size();
        v.insert(v.end(), L.v.begin(), L.v.end());
        s.reserve(s.size() + L.s.size());
        for ( size_t i = 0; i < L.s.size(); ++i )
            s.push_back(Chain(v0 + L.s[i].start, L.s[i].loop));
        
        // update bounds
        ext.add(L.ext);
        
        return *this;
    }
    
    // getters
    
    // get number of vertices.
    size_t size () const { return v.size(); }
    
    // get number of chains.
    size_t nchains () const { return s.size(); }
    
    // get the chain index for the chain that contains vertex i;
    // if no chain contains vertex i, returns the number of chains.
    size_t findChain ( size_t i ) const
    {
        if ( v.size() == 0 ) return 0;
        
        // define pivots
        size_t jp = s.size();
        size_t b = 0;
        size_t u = s.size();
        size_t j = u / 2;
        
        // find by bisection
        while ( j != jp ) {
            if ( s[j].start < i ) {
                if ( j+1 == s.size() ) return s.size();
                else b = j+1;
            }
            else if ( i < s[j].start ) u = j;
            else return j;
            jp = j;
            j = (b + u) / 2;
        }
        
        return j;
    }
    
    // get starting index of chain i.
    size_t chain ( size_t i ) const { return s[i].start; }
    
    // get whether chain i loops or not.
    bool loop ( size_t i ) const { return s[i].loop; }
    
    // get the object's extent.
    const Extent& extent () const { return ext; }
    
    // get the object's vertex vector.
    const std::vector<Point2d>& vertices () const { return v; }
    
    // basic operations
    
    // specify that a chain begins at vertex i and loops if c is
    // true; returns the chain number. If a chain already begins
    // at i, overwrites it.
    size_t addChainBreak ( size_t i, bool c )
    {
        size_t j = findChain(i);
        if ( j < s.size() && s[j].start == i ) s[j].loop = c;
        else s.insert(s.begin() + j, Chain(i, c));
        return j;
    }
    
    // remove the i'th break between two chains. If either was
    // looping, they are no longer considered looping
    // unless the beginning and end of the new chain are the same.
    void removeChainBreak ( size_t i )
    {
        if ( i <= 0 || i >= s.size() ) return;
        size_t j = i-1;
        size_t k = (i < s.size() ? s[i+1].start-1 : v.size()-1);
        
        // remove loop unless the start of j and end of i are same
        if ( s[j].loop || s[i].loop )
            s[j].loop = v[s[j].start] == v[k-1];
        
        // remove chain break
        s.erase(s.begin() + i);
    }
    
    // add a new chain before chain i that loops if c is true; let i
    // be greater than or equal to the number of chains to append
    // at the end of the polyline.
    void addChain ( size_t i, const std::vector<Point2d>& L, bool c )
    {
        if ( L.size() == 0 ) return;
        
        size_t v0 = v.size();
        size_t j0 = v0;
        if ( i >= s.size() ) i = s.size();
        else j0 = s[i].start;
        
        v.insert(v.begin() + j0, L.begin(), L.end()); // add vertices
        s.insert(s.begin() + i, Chain(j0, c)); // add chain break
        
        // update chain breaks
        ++i;
        for ( ; i < s.size(); ++i ) s[i].start += L.size();
        
        // update bounds
        if ( v0 == 0 ) ext = Extent(L[0].x, L[0].y, L[0].x, L[0].y);
        for ( size_t j = 0; j < L.size(); ++j )
            ext.add(Extent(L[j].x, L[j].y, L[j].x, L[j].y));
    }
    
    // create a polyline of chain i
    Polyline2d extractChain ( size_t i ) const
    {
        if ( i >= s.size() ) return Polyline2d();
        
        size_t j0 = s[i].start;
        size_t j1 = (i+1 < s.size() ? s[i+1].start : v.size());
        Polyline2d out;
        out.v.insert(out.v.end(), v.begin() + j0, v.begin() + j1);
        out.s.push_back(Chain(0, s[i].loop));
        out.computeBounds();
        
        return out;
    }
    
    // delete chain i from polyline
    void removeChain ( size_t i )
    {
        if ( i >= s.size() ) return;
        
        // remove vertices
        size_t j0 = s[i].start;
        size_t j1 = (i+1 < s.size() ? s[i+1].start : v.size());
        v.erase(v.begin() + j0, v.begin() + j1);
        
        // remove chain break
        s.erase(s.begin() + i);
        
        // update remaining chain breaks
        for ( ; i < s.size(); ++i ) s[i].start -= (j1 - j0);
        
        // recompute extent
        computeBounds();
    }
    
    // append to the head of a chain. If i is not a valid chain,
    // creates a new chain instead.
    void append_front ( size_t i, const std::vector<Point2d>& L,
        bool c )
    {
        if ( L.size() == 0 ) return;
        addChain(i, L, c);
        removeChainBreak(i+1);
    }
    
    // append to the tail of a chain. If i is not a valid chain,
    // adds a new chain instead.
    void append_back ( size_t i, const std::vector<Point2d>& L,
        bool c )
    {
        if ( L.size() == 0 ) return;
        addChain(i+1, L, c);
        removeChainBreak(i+1);
    }
    
    // compute minimum and maximum bounding points; ensure chains
    // are represented in data structure
    void computeBounds ()
    {
        if ( v.size() > 0 ) {
            if ( s.size() == 0 ) s.push_back(Chain(0, false));
            
            ext = Extent(v[0].x, v[0].y, v[0].x, v[0].y);
            for ( size_t i = 1; i < v.size(); ++i )
                ext.add(Extent(v[i].x, v[i].y, v[i].x, v[i].y));
        }
        else {
            s.clear();
            ext = Extent(2);
        }
    }
}; // Polyline2d


#endif // YOUNG_GEOMETRY_POLYLINE_20210111