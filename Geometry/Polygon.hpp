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
/////////////////////////////////////////////////////////////////////

#ifndef YOUNG_GEOMETRY_POLYGON_20201229
#define YOUNG_GEOMETRY_POLYGON_20201229

#include <vector>   // std::vector
#include "Line.hpp" // Point2d, Line2d
#include "../BYstdlib/sort.h" // quickSort


class Polygon2d {
public:
    std::vector<Point2d> v; // vertices
    Point2d bMin, bMax;     // bounds
    
    // constructors, destructor
    Polygon2d (
        const std::vector<Point2d>& V=std::vector<Point2d>() ) : v(V)
    {
        computeBounds();
    }
    Polygon2d ( const Polygon2d& G ) : v(G.v), bMin(G.bMin),
        bMax(G.bMax) {}
    ~Polygon2d () {}
    
    // operators
    Polygon2d& operator= ( const Polygon2d& G )
    {
        if ( &G == this ) return *this;
        v = G.v;
        bMin = G.bMin;
        bMax = G.bMax;
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
    bool operator< ( const Polygon2d& G ) const
    {
        if ( bMin < G.bMin ) return true;
        if ( bMin > G.bMin ) return false;
        if ( bMax < G.bMax ) return true;
        if ( bMax > G.bMax ) return false;
        return v.size() < G.v.size();
    }
    bool operator== ( const Polygon2d& G ) const
    {
        size_t i = 0;
        for ( ; i < v.size() && i < G.v.size(); ++i )
            if ( v[i] != G.v[i] ) return false;
        return !(i < v.size() || i < G.v.size());
    }
    bool operator!= ( const Polygon2d& G ) const
    {
        size_t i = 0;
        for ( ; i < v.size() && i < G.v.size(); ++i )
            if ( v[i] != G.v[i] ) return true;
        return i < v.size() || i < G.v.size();
    }
    
    // getters
    size_t size () const
    {
        return v.size();
    }
    
    // basic operations
    void computeBounds ()
    {
        if ( v.size() > 0 ) {
            bMin = bMax = v[0];
            for ( size_t i = 1; i < v.size(); ++i ) {
                if ( v[i].x < bMin.x ) bMin.x = v[i].x;
                else if ( v[i].x > bMax.x ) bMax.x = v[i].x;
                if ( v[i].y < bMin.y ) bMin.y = v[i].y;
                else if ( v[i].y > bMax.y ) bMax.y = v[i].y;
            }
        }
        else bMin = bMax = Point2d();
    }
}; // Polygon


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
//        or hole must be the same.
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
    double xmin, xmax, ymin, ymax;
    if ( s.a.x < s.b.x ) {
        xmin = s.a.x;
        xmax = s.b.x;
    }
    else {
        xmin = s.b.x;
        xmax = s.a.x;
    }
    if ( s.a.y < s.b.y ) {
        ymin = s.a.y;
        ymax = s.b.y;
    }
    else {
        ymin = s.b.y;
        ymax = s.a.y;
    }
    
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
        if ( (poly[i].x > xmax && poly[i+1].x > xmax)
                || (poly[i].x < xmin && poly[i+1].x < xmin)
                || (poly[i].y > ymax && poly[i+1].y > ymax)
                || (poly[i].y < ymin && poly[i+1].y < ymin) )
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
    
    // measure distance to each line segment
    size_t I, j = 0;
    for ( size_t i = 0; i+1 < poly.size(); ++i ) {
        // start of polygon
        if ( i == j ) I = i;
        
        // end of polygon
        else if ( poly[i] == poly[I] ) {
            j = i+1;
            continue;
        }
        
        // measure distance to line segment
        double sqd = sqDistToLine2d(p, Line2d(poly[i], poly[i+1]));
        if ( sqdmin < 0.0 || sqd < sqdmin ) sqdmin = sqd;
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
    
    // measure distance to each line segment
    size_t I, j = 0;
    for ( size_t i = 0; i+1 < poly.size(); ++i ) {
        // start of polygon
        if ( i == j ) I = i;
        
        // end of polygon
        else if ( poly[i] == poly[I] ) {
            j = i+1;
            continue;
        }
        
        // check if line segments intersect
        double sqd = sqDistToLine(s, Line2d(poly[i], poly[i+1]));
        if ( sqd < 0.0000001 ) return 0.0;
        if ( sqdmin < 0.0 || sqd < sqdmin ) sqdmin = sqd;
    }
    
    return sqdmin;
}

double distToEdge2d ( const Polygon2d& poly, const Line2d& s )
{
    return sqrt(sqDistToEdge2d(poly, s));
}


#endif // YOUNG_GEOMETRY_POLYGON_20201229