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
// 02/05/2021 - Brennan Young                                      //
// - reorganized to define methods outside of class definition.    //
// - added chainEnd.                                               //
// 02/08/2021 - Brennan Young                                      //
// - corrected error in findChain.                                 //
// 02/12/2021 - Brennan Young                                      //
// - findChain now correctly updates the upper bound as it         //
//   iterates.                                                     //
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
    // constructors, destructor
    Polyline2d(const std::vector<Point2d>&, bool);
    Polyline2d(const Polygon2d&);
    Polyline2d(const Polyline2d&);
    ~Polyline2d();
    
    // operators
    Polyline2d& operator=(const Polyline2d&);   // assignment
    Point2d& operator[](size_t);                // element access
    const Point2d& operator[](size_t) const;    // element access
    bool operator<(const Polyline2d&) const;    // less-than
    bool operator==(const Polyline2d&) const;   // equal-to
    bool operator!=(const Polyline2d&) const;   // not equal-to
    Polyline2d operator+(const Polyline2d&) const; // append
    Polyline2d& operator+=(const Polyline2d&);  // append
    
    // getters
    size_t size() const;                    // number of vertices
    size_t nchains() const;                 // number of chains
    size_t findChain(size_t) const;         // chain containing i
    size_t chain(size_t) const;             // 1st vertex 0 in chain
    size_t chainEnd(size_t) const;          // last+1 vertex in chain
    bool loop(size_t) const;                // if chain loops
    const Extent& extent() const;           // get extent object
    const std::vector<Point2d>& vertices() const; // vertex vector
    
    // basic operations
    size_t addChainBreak(size_t, bool);     // new chain break
    void removeChainBreak(size_t);          // remove chain break
    void addChain(size_t, const std::vector<Point2d>&, bool);
    Polyline2d extractChain(size_t) const;  // get chain vertices
    void removeChain(size_t);               // delete chain vertices
    void insert(size_t, const Point2d&);    // add vertex
    template <class Iterator> void insert(  // add vertices
        size_t, const Iterator&, const Iterator&);
    void push_back(const Point2d&);         // append vertex to tail
    void erase(size_t);                     // delete vertex
    
    // append to the head/tail of a chain
    void append_front(size_t, const std::vector<Point2d>&, bool);
    void append_back(size_t, const std::vector<Point2d>&, bool);
    
    void computeBounds();                   // determine extent
}; // Polyline2d


// CONSTRUCTORS / DESTRUCTOR ////////////////////////////////////////


// constructor
Polyline2d::Polyline2d (
    const std::vector<Point2d>& V=std::vector<Point2d>(),
    bool c=false )
: v(V), ext(Extent(2))
{
    s.push_back(Chain(0, c));
    computeBounds();
}

// construct from polygon
Polyline2d::Polyline2d ( const Polygon2d& G )
: v(G.vertices()), ext(G.extent())
{
    for ( size_t i = 0; i < G.nchains(); ++i )
        s.push_back(Chain(G.chain(i), G.loop(i)));
}

// copy constructor
Polyline2d::Polyline2d ( const Polyline2d& L )
: v(L.v), s(L.s), ext(L.ext)
{}

// destructor
Polyline2d::~Polyline2d () {}


// OPERATORS ////////////////////////////////////////////////////////


// assignment
Polyline2d& Polyline2d::operator= ( const Polyline2d& L )
{
    if ( &L == this ) return *this;
    v = L.v;
    s = L.s;
    ext = L.ext;
    return *this;
}

// element access
Point2d& Polyline2d::operator[] ( size_t i )
{
    return v[i];
}

const Point2d& Polyline2d::operator[] ( size_t i ) const
{
    return v[i];
}

// less-than comparison
bool Polyline2d::operator< ( const Polyline2d& L ) const
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

// equal-to comparison
bool Polyline2d::operator== ( const Polyline2d& L ) const
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

// not-equal-to comparison
bool Polyline2d::operator!= ( const Polyline2d& L ) const
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
Polyline2d Polyline2d::operator+ ( const Polyline2d& L ) const
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
Polyline2d& Polyline2d::operator+= ( const Polyline2d& L )
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


// GETTERS //////////////////////////////////////////////////////////


// get number of vertices.
size_t Polyline2d::size () const { return v.size(); }

// get number of chains.
size_t Polyline2d::nchains () const { return s.size(); }

// get the chain index for the chain that contains vertex i;
// if no chain contains vertex i, returns the number of chains
size_t Polyline2d::findChain ( size_t i ) const
{
    if ( v.size() == 0 ) return 0;
    
    // define pivots
    size_t jp = s.size();
    size_t b = 0;        // lower pivot
    size_t u = s.size(); // upper pivot
    size_t j = u / 2;    // chain
    size_t i1 = j+1 < s.size() ? s[j+1].start : v.size(); // up bound
    
    // find by bisection
    while ( j != jp ) {
        if ( i >= i1 ) {
            if ( j+1 == s.size() ) return s.size();
            else b = j+1;
        }
        else if ( i < s[j].start ) u = j;
        else return j;
        jp = j;
        j = (b + u) / 2;
        i1 = j+1 < s.size() ? s[j+1].start : v.size(); // new bound
    }
    
    return j;
}

// get starting index of chain i
size_t Polyline2d::chain ( size_t i ) const { return s[i].start; }

// get ending index of chain i (not inclusive)
size_t Polyline2d::chainEnd ( size_t i ) const
{
    return i+1 < s.size() ? s[i+1].start : v.size();
}

// get whether chain i loops or not
bool Polyline2d::loop ( size_t i ) const { return s[i].loop; }

// get the object's extent
const Extent& Polyline2d::extent () const { return ext; }

// get the object's vertex vector
const std::vector<Point2d>& Polyline2d::vertices () const
{
    return v;
}


// BASIC OPERATIONS /////////////////////////////////////////////////


// specify that a chain begins at vertex i and loops if c is
// true; returns the chain number. If a chain already begins
// at i, overwrites it.
size_t Polyline2d::addChainBreak ( size_t i, bool c )
{
    size_t j = findChain(i);
    if ( j < s.size() && s[j].start == i ) s[j].loop = c;
    else s.insert(s.begin() + j, Chain(i, c));
    return j;
}

// remove the i'th break between two chains. If either was
// looping, they are no longer considered looping
// unless the beginning and end of the new chain are the same.
void Polyline2d::removeChainBreak ( size_t i )
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
// at the end of the object.
void Polyline2d::addChain ( size_t i, const std::vector<Point2d>& L,
    bool c )
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
Polyline2d Polyline2d::extractChain ( size_t i ) const
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

// delete chain i from object
void Polyline2d::removeChain ( size_t i )
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

// insert a vertex
void Polyline2d::insert ( size_t i, const Point2d& p )
{
    size_t v0 = v.size();
    size_t j = findChain(i);
    v.insert(v.begin() + i, p);
    ++j;
    for ( ; j < s.size(); ++j ) s[j].start += 1;
    if ( v0 == 0 ) ext = Extent(p.x, p.y, p.x, p.y);
    else ext.add(Extent(p.x, p.y, p.x, p.y));
}

// insert vertices
template <class Iterator>
void Polyline2d::insert ( size_t i, const Iterator& begin,
    const Iterator& end )
{
    size_t v0 = v.size();
    v.insert(v.begin() + i, begin, end);
    Iterator it = begin;
    if ( v0 == 0 ) {
        ext = Extent(it->x, it->y, it->x, it->y);
        ++it;
    }
    for ( ; it != end; ++it )
        ext.add(Extent(it->x, it->y, it->x, it->y));
}

// append vertex to tail
void Polyline2d::push_back ( const Point2d& p )
{
    if ( s.size() == 0 ) s.push_back(Chain(0,false));
    
    v.push_back(p);
    if ( v.size() == 1 ) ext = Extent(p.x, p.y, p.x, p.y);
    else ext.add(Extent(p.x, p.y, p.x, p.y));
}

// delete vertex
void Polyline2d::erase ( size_t i )
{
    size_t j = findChain(i);
    size_t j0 = j;
    v.erase(v.begin() + i);
    ++j;
    for ( ; j < s.size(); ++j ) s[j].start -= 1;
    if ( j0+1 > s.size() && s[j0].start == s[j0+1].start )
        s.erase(s.begin() + j0);
    computeBounds();
}

// append to the head of a chain. If i is not a valid chain,
// creates a new chain instead.
void Polyline2d::append_front ( size_t i,
    const std::vector<Point2d>& L, bool c )
{
    if ( L.size() == 0 ) return;
    addChain(i, L, c);
    removeChainBreak(i+1);
}

// append to the tail of a chain. If i is not a valid chain,
// adds a new chain instead.
void Polyline2d::append_back ( size_t i,
    const std::vector<Point2d>& L, bool c )
{
    if ( L.size() == 0 ) return;
    addChain(i+1, L, c);
    removeChainBreak(i+1);
}

// compute minimum and maximum bounding points; ensure chains
// are represented in data structure
void Polyline2d::computeBounds ()
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


#endif // YOUNG_GEOMETRY_POLYLINE_20210111