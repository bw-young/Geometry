/////////////////////////////////////////////////////////////////////
// Objects and methods for managing geometric polygons (2d).       //
/////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////
// -- HISTORY ---------------------------------------------------- //
// 12/29/2020 - Brennan Young                                      //
// - created                                                       //
// - migrated inPoly from ObjectPolyOps.h.                         //
// 01/07/2020 - Brennan Young                                      //
// - migrated segmentPoly_intersect2d from ObjectPolyOps.h.        //
// - migrated sqDistToEdge and distToEdge from ObjectPolyOps.h.    //
// 01/22/2021 - Brennan Young                                      //
// - replaced Point2d bMin and bMax objects with Extent ext.       //
// 02/05/2021 - Brennan Young                                      //
// - restructure to contain information about where each vertex    //
//   chain (sub-polygon) begins and to protect object members.     //
// - reorganized to define methods outside of class definition.    //
// 02/08/2021 - Brennan Young                                      //
// - corrected error in findChain.                                 //
/////////////////////////////////////////////////////////////////////

#ifndef YOUNG_GEOMETRY_POLYGON_20201229
#define YOUNG_GEOMETRY_POLYGON_20201229

#include <vector>             // std::vector
#include "Extent.hpp"         // Extent
#include "Line.hpp"           // Point2d, Line2d
#include "../BYstdlib/sort.h" // quickSort


class Polygon2d {
private:
    std::vector<Point2d> v; // vertices
    std::vector<size_t> s;  // chain breaks
    Extent ext;             // bounds
    
public:
    // constructors, destructor
    Polygon2d(const std::vector<Point2d>&);
    Polygon2d(const Polygon2d&);
    ~Polygon2d();
    
    // operators
    Polygon2d& operator=(const Polygon2d&);     // assignment
    Point2d& operator[](size_t i);              // element access
    const Point2d& operator[]( size_t) const;   // element access
    bool operator<(const Polygon2d&) const;     // less-than
    bool operator==(const Polygon2d&) const;    // equal-to
    bool operator!=(const Polygon2d&) const;    // not-equal-to
    
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
    void addChain(size_t, const std::vector<Point2d>&); // new chain
    Polygon2d extractChain(size_t) const;   // get chain vertices
    void removeChain(size_t);               // delete chain vertices
    void insert(size_t, const Point2d&);    // add vertex
    template <class Iterator> void insert(  // add vertices
        size_t, const Iterator&, const Iterator&);
    void push_back(const Point2d&);         // append vertex to tail
    void erase(size_t);                     // delete vertex
    void computeBounds();                   // determine extent
}; // Polygon


// CONSTRUCTORS / DESTRUCTOR ////////////////////////////////////////


// constructor
Polygon2d::Polygon2d (
    const std::vector<Point2d>& V=std::vector<Point2d>() )
: v(V), ext(Extent(2))
{
    // identify closed chains
    size_t i0=0, i1=0;
    for ( size_t i = 0; i < v.size(); ++i ) {
        // start of polygon
        if ( i == i1 ) {
            i0 = i1;
            s.push_back(i);
        }
        
        // end of polygon
        else if ( i == i0 ) i1 = i+1;
    }
    
    computeBounds();
}

// copy constructor
Polygon2d::Polygon2d ( const Polygon2d& G )
: v(G.v), s(G.s), ext(G.ext)
{}

// destructor
Polygon2d::~Polygon2d () {}


// OPERATORS ////////////////////////////////////////////////////////


// assignment
Polygon2d& Polygon2d::operator= ( const Polygon2d& G )
{
    if ( &G == this ) return *this;
    v    = G.v;
    s    = G.s;
    ext  = G.ext;
    return *this;
}

// element access
Point2d& Polygon2d::operator[] ( size_t i )
{
    return v[i];
}

const Point2d& Polygon2d::operator[] ( size_t i ) const
{
    return v[i];
}

// less-than comparison
bool Polygon2d::operator< ( const Polygon2d& G ) const
{
    // sort by bounding box first
    if ( ext < G.ext ) return true;
    if ( ext > G.ext ) return false;
    
    // then sort by number of vertices
    if ( v.size() < G.v.size() ) return true;
    if ( v.size() > G.v.size() ) return false;
    
    // then sort by vertices
    for ( size_t i = 0; i < v.size(); ++i ) {
        if ( v[i] < G.v[i] ) return true;
        if ( v[i] > G.v[i] ) return false;
    }
    return false;
}

// equal-to comparison
bool Polygon2d::operator== ( const Polygon2d& G ) const
{
    size_t i = 0;
    for ( ; i < v.size() && i < G.v.size(); ++i )
        if ( v[i] != G.v[i] ) return false;
    return !(i < v.size() || i < G.v.size());
}

// not-equal-to comparison
bool Polygon2d::operator!= ( const Polygon2d& G ) const
{
    size_t i = 0;
    for ( ; i < v.size() && i < G.v.size(); ++i )
        if ( v[i] != G.v[i] ) return true;
    return i < v.size() || i < G.v.size();
}


// GETTERS //////////////////////////////////////////////////////////

// get number of vertices
size_t Polygon2d::size () const { return v.size(); }

// get number of chains
size_t Polygon2d::nchains () const { return s.size(); }

// get the chain index for the chain that contains vertex i;
// if no chain contains vertex i, returns the number of chains
size_t Polygon2d::findChain ( size_t i ) const
{
    if ( v.size() == 0 ) return 0;
    
    // define pivots
    size_t jp = s.size();
    size_t b = 0;        // lower pivot
    size_t u = s.size(); // upper pivot
    size_t j = u / 2;    // chain
    size_t i1 = j+1 < s.size() ? s[j+1] : v.size(); // upper bound
    
    // find by bisection
    while ( j != jp ) {
        if ( i >= i1 ) {
            if ( j+1 == s.size() ) return s.size();
            else b = j+1;
        }
        else if ( i < s[j] ) u = j;
        else return j;
        jp = j;
        j = (b + u) / 2; // new chain
        i1 = j+1 < s.size() ? s[j+1] : v.size(); // new upper bound
    }
    
    return j;
}

// get starting index of chain i
size_t Polygon2d::chain ( size_t i ) const { return s[i]; }

// get ending index of chain i (not inclusive)
size_t Polygon2d::chainEnd ( size_t i ) const
{
    return i+1 < s.size() ? s[i+1] : v.size();
}

// get whether chain i loops or not (always true for polygons)
bool Polygon2d::loop ( size_t i ) const { return true; }

// get the object's extent
const Extent& Polygon2d::extent () const { return ext; }

// get the object's vertex vector
const std::vector<Point2d>& Polygon2d::vertices () const
{
    return v;
}


// BASIC OPERATIONS /////////////////////////////////////////////////

// add a new chain to the object
// -- Arguments --
// i : the chain index to insert the new chain. If i >= number of
//     chains, appens a new chain at the end.
// G : vertices in the new chain. If G.back() != G[0], pushes G[0]
//     onto the tail of the new chain.
// -- Returns --
// Nothing.
void Polygon2d::addChain ( size_t i, const std::vector<Point2d>& G )
{
    if ( G.size() == 0 ) return;
    
    size_t v0 = v.size();
    size_t j0 = v0;
    if ( i >= s.size() ) i = s.size();
    else j0 = s[i];
    
    bool closed = G.back() == G[0];
    int Gsize = G.size() + closed ? 0 : 1;
    
    v.insert(v.begin() + j0, G.begin(), G.end()); // add vertices
    if ( !closed ) v.insert(v.begin() + j0 + G.size(), G[0]);
    s.insert(s.begin() + i, j0); // add chain break
    
    // update chain breaks
    ++i;
    for ( ; i < s.size(); ++i ) s[i] += Gsize;
    
    // update bounds
    if ( v0 == 0 ) ext = Extent(G[0].x, G[0].y, G[0].x, G[0].y);
    for ( size_t j = 0; j < G.size(); ++j )
        ext.add(Extent(G[j].x, G[j].y, G[j].x, G[j].y));
}

// create a polygon of chain i
Polygon2d Polygon2d::extractChain ( size_t i ) const
{
    if ( i >= s.size() ) return Polygon2d();
    
    size_t j0 = s[i];
    size_t j1 = (i+1 < s.size() ? s[i+1] : v.size());
    Polygon2d out;
    out.v.insert(out.v.end(), v.begin() + j0, v.begin() + j1);
    out.s.push_back(0);
    out.computeBounds();
    
    return out;
}

// delete chain i from object
void Polygon2d::removeChain ( size_t i )
{
    if ( i >= s.size() ) return;
    
    // remove vertices
    size_t j0 = s[i];
    size_t j1 = (i+1 < s.size() ? s[i+1] : v.size());
    v.erase(v.begin() + j0, v.begin() + j1);
    
    // remove chain break
    s.erase(s.begin() + i);
    
    // update remaining chain breaks
    for ( ; i < s.size(); ++i ) s[i] -= (j1 - j0);
    
    // recompute extent
    computeBounds();
}

// insert a vertex
void Polygon2d::insert ( size_t i, const Point2d& p )
{
    size_t v0 = v.size();
    size_t j = findChain(i);
    v.insert(v.begin() + i, p);
    ++j;
    for ( ; j < s.size(); ++j ) s[j] += 1;
    if ( v0 == 0 ) ext = Extent(p.x, p.y, p.x, p.y);
    else ext.add(Extent(p.x, p.y, p.x, p.y));
}

// insert vertices
template <class Iterator>
void Polygon2d::insert ( size_t i, const Iterator& begin,
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

// append vertex to tail; creates a new chain if adds on to a
// completed one
void Polygon2d::push_back ( const Point2d& p )
{
    // new chain?
    if ( s.size() == 0 ) s.push_back(0);
    else if ( v.size() - s.back() > 2 && v.back() == v[s.back()] )
        s.push_back(v.size());
    
    v.push_back(p);
    
    // update extent
    if ( v.size() == 1 ) ext = Extent(p.x, p.y, p.x, p.y);
    else ext.add(Extent(p.x, p.y, p.x, p.y));
}

// delete vertex
void Polygon2d::erase ( size_t i )
{
    size_t j = findChain(i);
    size_t j0 = j;
    v.erase(v.begin() + i);
    ++j;
    for ( ; j < s.size(); ++j ) s[j] -= 1;
    if ( j0+1 > s.size() && s[j0] == s[j0+1] )
        s.erase(s.begin() + j0);
    computeBounds();
}

// determine bounds
void Polygon2d::computeBounds ()
{
    if ( v.size() > 0 ) {
        if ( s.size() == 0 ) s.push_back(0);
        ext = Extent(v[0].x, v[0].y, v[0].x, v[0].y);
        for ( size_t i = 1; i < v.size(); ++i )
            ext.add(Extent(v[i].x, v[i].y, v[i].x, v[i].y));
    }
    else ext = Extent(2);
}



/////////////////////////////////////////////////////////////////////
// Primitive Operations                                            //
/////////////////////////////////////////////////////////////////////


// Check if given location (x,y) is inside polygon using the winding
// number method.
// -- Source --
// - Dan Sunday, 2001, Inclusion of a point in a polygon:,
//   http://geomalgorithms.com/a03-_inclusion.html (accessed
//   4/16/2020)
// -- Arguments --
// poly : polygon to check. The first and last point in each polygon
//        or hole (vertex chain) must be the same.
// pt   : point representing location to check. Point object must
//        have public real or integer x and y members.
// -- Returns --
// The winding number, or the number of times that the polygon wraps
// around (x,y). Winding number is > 0 if (x,y) is inside the
// polygon. For a polygon with crossing edges, an even winding number
// indicates that (x,y) is inside one of the "internal" parts of
// the polygon, which would be identified as a hole with the ray
// casting algorithm. Therefore:
//   int w = inPoly(poly,pt);
//   bool isInPolygon = w > 0 && w % 2 != 0;
template <class Point>
int inPoly ( const Polygon2d& poly, Point pt )
{
    if ( poly.size() == 0 ) return 0;
    
    double left;
    size_t I, i, j = 0;
    int n = 0;
    for ( i = 0; i < poly.size() - 1; ++i ) {
        // start of polygon
        if ( i == j ) I = i;
        
        // end of polygon
        else if ( poly[i] == poly[I] ) {
            j = i+1;
            continue;
        }
        
        // check if the segment crosses the +x ray
        if ( poly[i].y <= pt.y ) {
            // upward crossing
            if ( poly[i+1].y > pt.y ) {
                left = (poly[i+1].x - poly[i].x)
                    * (pt.y - poly[i].y)
                    - (pt.x - poly[i].x)
                    * (poly[i+1].y - poly[i].y);
                if ( left > 0.0 ) ++n;
            }
        }
        else {
            // downward crossing
            if ( poly[i+1].y <= pt.y ) {
                left = (poly[i+1].x - poly[i].x)
                    * (pt.y - poly[i].y)
                    - (pt.x - poly[i].x)
                    * (poly[i+1].y - poly[i].y);
                if ( left < 0.0 ) --n;
            }
        }
    }
    
    return n;
}

// Determine if and where a ray intersects a polygon.
// -- Arguments --
// r0   : ray origin.
// r1   : point representing the ray direction vector.
// poly : polygon object.
// rmax : maximum distance from r0 to include.
// -- Returns --
// A vector of the distances from the ray origin to the points of
// intersection, in order of increasing distance. This is empty if
// they do not intersect.
template <class Point>
std::vector<double> rayPoly_intersect2d ( const Point& r0,
    const Point& r1, const Polygon2d& poly, double rmax )
{
    std::vector<double> intersections;
    
    double d, dx, dy, dot;
    bool d0, dd, R, R0, Rp;
    d0 = dd = R = R0 = Rp = false;
    
    // scan for edge of object
    size_t m, I, j = 0;
    for ( size_t i = 0; i < poly.size() - 1; ++i ) {
        // start of polygon
        if ( i == j ) {
            I = j;
            m = -1;
        }
        
        // end of polygon
        else if ( poly[i] == poly[I] ) {
            j = i+1;
            continue;
        }
        
        // determine if segment's bounding box intersects ray's
        if ( (r1.x < 0.0 && poly[i].x > r0.x && poly[i+1].x > r0.x)
          || (r1.x > 0.0 && poly[i].x < r0.x && poly[i+1].x < r0.x)
          || (r1.y < 0.0 && poly[i].y > r0.y && poly[i+1].y > r0.y)
          || (r1.y > 0.0 && poly[i].y < r0.y && poly[i+1].y < r0.y) )
            continue;
        
        // get distance to point of intersection
        d = raySegment_intersect2d(r0, r1, poly[i], poly[i+1]);
        
        // skip if doesn't intersect segment
        if ( d < 0.0 ) continue;
        
        // determine if same as first or previous intersection
        d0 = poly[i+1] == poly[I]
            && m < intersections.size()
            && fabs(intersections[m] - d) < 0.0000001;
        dd = intersections.size() > 0
            && fabs(intersections.back() - d) < 0.0000001;
        
        // get whether segment direction is clockwise of ray
        dx = poly[i+1].x - poly[i].x; // component vector of segment
        dy = poly[i+1].y - poly[i].y;
        
        // get dot product of right-orthogonal vector (+ = clockwise)
        dot = r1.x * (-dy) + r1.y * dx;
        R = dot > 0.0;
        
        if ( i == I ) R0 = Rp = R;
        
        // record intersection only if not duplicating a vertex
        // unless ray touches but doesn't cross polygon boundary
        if ( !((d0 && R0 == R) || (dd && Rp == R)) ) {
            intersections.push_back(d);
            if ( i == I ) {
                m = intersections.size() - 1;
                R0 = R;
            }
        }
        
        Rp = R;
    }
    
    bystd::quickSort(&intersections);
    return intersections;
}

// Determine if and where a line segment intersects a polygon.
// -- Arguments --
// poly : the polygon being evaluated.
// s    : the line segment being evaluated.
// -- Returns --
// A sorted vector of intersection distances from s0 toward s1.
std::vector<double> segmentPoly_intersect2d (
    const Polygon2d& poly, const Line2d& s )
{
    std::vector<double> intersections;
    if ( poly.size() == 0 ) return intersections;
    
    // segment extent
    Extent sExt (s.a.x, s.a.y, s.b.x, s.b.y, 1.0);
    if ( !poly.extent().overlap(sExt) ) return intersections;
    
    double dot;
    bool d0, dd, R, R0, Rp;
    d0 = dd = R = R0 = Rp = false;
    
    // scan for edge of object
    double d, dx, dy;
    size_t m, I = 0, j = 0;
    for ( size_t i = 0; i < poly.size() - 1; ++i ) {
        // start of polygon
        if ( i == j ) {
            I = j;
            m = -1;
        }
        
        // end of polygon
        else if ( poly[i] == poly[I] ) {
            j = i+1;
            continue;
        }
        
        // determine if segments' bounding boxes don't intersect
        if ( (poly[i].x > sExt.xmin()
                    && poly[i+1].x > sExt.xmax())
                || (poly[i].x < sExt.xmin()
                    && poly[i+1].x < sExt.xmin())
                || (poly[i].y > sExt.ymax()
                    && poly[i+1].y > sExt.ymax())
                || (poly[i].y < sExt.ymin()
                    && poly[i+1].y < sExt.ymin()) )
            continue;
        
        // get distance to point of intersection
        d = raySegment_intersect2d(s.a, s.r, poly[i], poly[i+1]);
        
        // skip if doesn't intersect segment
        if ( d < 0.0 || d > s.len ) continue;
        
        // determine if same as first or previous intersection
        d0 = poly[i+1] == poly[I]
            && m < intersections.size()
            && fabs(intersections[m] - d) < 0.0000001;
        dd = intersections.size() > 0
            && fabs(intersections.back() - d) < 0.0000001;
        
        // get whether segment direction is clockwise of ray
        dx = poly[i+1].x - poly[i].x; // component vector of segment
        dy = poly[i+1].y - poly[i].y;
        
        // get dot product of right-orthogonal vector (- = clockwise)
        dot = s.r.x * (-dy) + s.r.y * dx;
        R = dot > 0.0;
        
        if ( i == I ) R0 = Rp = R;
        
        // record intersection only if not duplicating a vertex
        // unless ray touches but doesn't cross polygon boundary
        if ( !((d0 && R0 == R) || (dd && Rp == R)) ) {
            intersections.push_back(d);
            if ( i == I ) {
                m = intersections.size() - 1;
                R0 = R;
            }
        }
        
        Rp = R;
    }
    
    bystd::quickSort(&intersections);
    return intersections;
}

// Return the distance between the given point and the nearest edge
// of the polygon
// -- Arguments --
// poly : Polygon object to evaluate.
// p    : Point object representing the location to check.
// -- Returns --
// The distance to the nearest polygon edge. Returns a negative
// number if no distance could be determined.
double sqDistToEdge2d ( const Polygon2d& poly, const Point2d& p )
{
    // get distance to first point
    double sqdmin = -1.0;
    
    // for each chain in object
    for ( size_t ii = 0; ii < poly.nchains(); ++ii ) {
        // measure distance to each line segment in chain
        size_t i1 = poly.chainEnd(ii);
        for ( size_t i = poly.chain(ii); i < i1; ++i ) {
            double sqd = sqDistToLine(p, Line2d(poly[i], poly[i+1]));
            if ( sqdmin < 0.0 || sqd < sqdmin ) sqdmin = sqd;
        }
    }
    
    return sqdmin;
}

double distToEdge2d ( const Polygon2d& poly, const Point2d& p )
{
    return sqrt(sqDistToEdge2d(poly, p));
}

// Return the distance between the given line segment and the nearest
// edge of the polygon.
// -- Arguments --
// poly : Polygon object to evaluate.
// s    : Line segment to evaluate.
// -- Returns --
// The distance to the nearest polygon edge. Returns a negative
// number if no distance could be determined.
double sqDistToEdge2d ( const Polygon2d& poly, const Line2d& s )
{
    double sqdmin = -1.0;
    
    // for each chain in object
    for ( size_t ii = 0; ii+1 < poly.nchains(); ++ii ) {
        // measure distance to each line segment in chain
        size_t i1 = poly.chainEnd(ii);
        for ( size_t i = poly.chain(ii); i+1 < i1; ++i ) {
            double sqd = sqDistToLine(s, Line2d(poly[i], poly[i+1]));
            if ( sqd < 0.0000001 ) return 0.0; // intersect
            if ( sqdmin < 0.0 || sqd < sqdmin ) sqdmin = sqd;
        }
    }
    
    return sqdmin;
}

double distToEdge2d ( const Polygon2d& poly, const Line2d& s )
{
    return sqrt(sqDistToEdge2d(poly, s));
}


#endif // YOUNG_GEOMETRY_POLYGON_20201229