/////////////////////////////////////////////////////////////////////
// Objects and methods for managing geometric lines.               //
/////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////
// -- HISTORY ---------------------------------------------------- //
// 04/23/2020 - Brennan Young                                      //
// - created                                                       //
// - migrated Point2d, Point3d, raySegment_intersect2d,            //
//   and rayPoly_intersect2d from GridObjPolyOps.h.                //
// 08/21/2020 - Brennan Young                                      //
// - added distToEdge.                                             //
// 11/02/2020 - Brennan Young                                      //
// - rayPoly_intersect2d and segmentPoly_intersect2d no longer     //
//   double-count intersected vertices.                            //
// 11/11/2020                                                      //
// - rayPoly_intersect2d now double-counts vertices where the ray  //
//   both enters and exits the polygon at that vertex.             //
// 11/18/2020                                                      //
// - rayPoly_intersect2d and segmentPoly_intersect2d vertex-       //
//   duplicate avoidance has been corrected.                       //
// 12/11/2020                                                      //
// - migrated Line2d and raySegment_intersect2d from               //
//   ObjectPolyOps.h.                                              //
// - added Line2d class.                                           //
// - added sqDistToLine2d, sqDistToLine, distToLine2d, and         //
//   distToLine.                                                   //
// 01/26/2021                                                      //
// - added angle2d.                                                //
// 02/04/2021                                                      //
// - moved angle2d to VectorOps.hpp.                               //
// - added lineIntersect2d.                                        //
// 02/05/2021                                                      //
// - lineIntersect2d now returns -1 intead of 0 when no            //
//   intersection.                                                 //
// 02/08/2021                                                      //
// - added sqMaxDistToLine and maxDistToLine functions.            //
// 02/09/2021                                                      //
// - migrate angle(Line) from VectorOps.hpp                        //
/////////////////////////////////////////////////////////////////////

#ifndef YOUNG_GEOMETRY_LINE_20201211
#define YOUNG_GEOMETRY_LINE_20201211

#include "VectorOps.hpp" // includes Point.hpp


// Basic 2D line object.
class Line2d {
public:
    Point2d a, b, r;    // start, end, and direction vector a->b
    double dx, dy, len; // euclidean distance from a to b
    
    // constructors, destructor
    Line2d ( const Point2d& A=Point2d(),
        const Point2d& B=Point2d() ) : a(A), b(B), dx(B.x-A.x),
        dy(B.y-A.y)
    {
        len = sqrt(dx*dx + dy*dy);
        if ( len > 0.0 ) r = Point2d(dx/len, dy/len);
    }
    Line2d ( const Line2d& L ) : a(L.a), b(L.b), r(L.r), dx(L.dx),
        dy(L.dy), len(L.len) {}
    ~Line2d () {}
    
    // operators
    Line2d& operator= ( const Line2d& L )
    {
        if ( &L == this ) return *this;
        a = L.a;
        b = L.b;
        r = L.r;
        dx = L.dx;
        dy = L.dy;
        len = L.len;
        return *this;
    }
    bool operator< ( const Line2d& L ) const
    {
        return (a < L.a && a < L.b) || (b < L.a && b < L.b);
    }
    bool operator== ( const Line2d& L ) const
    {
        return a == L.a && b == L.b;
    }
    bool operator!= ( const Line2d& L ) const
    {
        return a != L.a || b != L.b;
    }
};


/////////////////////////////////////////////////////////////////////
// Primitive Operations                                            //
/////////////////////////////////////////////////////////////////////


// Determine if and where a ray intersects a line segment.
// -- Arguments --
// r0 : ray origin.
// r1 : point representing the ray direction vector.
// p0 : line segment end-point.
// p1 : line segment end-point.
// -- Returns --
// The distance from the ray origin to the point of the intersection.
// This is negative if they do not intersect.
template <class Point>
double raySegment_intersect2d ( const Point& r0, const Point& r1,
    const Point& p0, const Point& p1 )
{
    // check for segment parallel to ray
    Point seg = Point(p1.x - p0.x, p1.y - p0.y); // direction vector
    Point segNrm = Point(seg.y, -seg.x); // vector normal to segment
    double d = r1.x * segNrm.x + r1.y * segNrm.y;
    if ( fabs(d) < 0.0000001 ) return -1.0;
    
    // identify where an intersection would occur
    // (t=distance from ray origin, s=distance from p0 to p1)
    seg = Point(p0.x - r0.x, p0.y - r0.y);
    double t = (seg.x * segNrm.x + seg.y * segNrm.y) / d;
    double s = (r1.y * seg.x - r1.x * seg.y) / d;
    
    if ( t >= 0.0 && s >= 0.0 && s <= 1.0 ) return t;
    return -1.0;
}

// Determine if and where two line segments intersect.
// -- Arguments --
// a : Line2d object to check for intersection.
// b : Line2d object to check for intersection.
// -- Returns --
// The distance from a.a toward a.b where a intersects b, or negative
// if they don't intersect.
double lineIntersect2d ( const Line2d& a, const Line2d& b )
{
    double d = raySegment_intersect2d(a.a, a.r, b.a, b.b);
    if ( d > a.len + 0.0000001 ) return -d;
    return d;
}

// Measure the distance from a point to a line segment.
// -- Arguments --
// p : point.
// s : line segment.
// -- Returns --
// The minimum distance from the point to the line segment.
template <class Point, class Line>
double sqDistToLine2d ( const Point& p, const Line& s )
{
    if ( s.len < 0.0000001 ) {
        // line segment of zero length - treat as a point
        return sqDist2d(p, s.a);
    }
    
    double u = ((p.x - s.a.x) * s.dx + (p.y - s.a.y) * s.dy) / s.len;
    
    if ( u > 1 ) u = 1;
    else if ( u < 0 ) u = 0;
    
    double xx = s.a.x + u * s.dx;
    double yy = s.a.y + u * s.dy;
    
    double dx = xx - p.x;
    double dy = yy - p.y;
    
    return dx*dx + dy*dy;
}
sqDistToLine ( const Point2d& p, const Line2d& s )
{
    return sqDistToLine2d(p, s);
}

template <class Point, class Line>
double distToLine2d ( const Point& p, const Line& s )
{
    return sqrt(sqDistToLine(p, s));
}
distToLine ( const Point2d& p, const Line2d& s )
{
    return distToLine2d(p, s);
}

// Measure the minimum distance between line segments.
// -- Arguments --
// a : line segment to compare.
// b : line segment to compare.
// -- Returns --
// The minimum square distance between the line segments.
template <class Line>
double sqDistToLine2d ( const Line& a, const Line& b )
{
    // if zero-length segment, treat as a point
    if ( a.len < 0.000001 ) return sqDistToLine2d(a.a, b);
    
    // check if line segments intersect
    double d = raySegment_intersect2d(a.a, a.r, b.a, b.b);
    if ( d >= 0.0 && d < a.len ) return 0.0;
    
    // check shortest sq distances between segment ends
    d = fabs(d);
    
    double dx = b.a.x - a.a.x;
    double dy = b.a.y - a.a.y;
    double len = dx*dx + dy*dy;
    if ( len < d ) d = len;
    
    dx = b.b.x - a.a.x;
    dy = b.b.y - a.a.y;
    len = dx*dx + dy*dy;
    if ( len < d ) d = len;
    
    dx = b.a.x - a.b.x;
    dy = b.a.y - a.b.y;
    len = dx*dx + dy*dy;
    if ( len < d ) d = len;
    
    dx = b.b.x - a.b.x;
    dy = b.b.y - a.b.y;
    len = dx*dx + dy*dy;
    if ( len < d ) d = len;
    
    return d;
}
double sqDistToLine ( const Line2d& a, const Line2d& b )
{
    return sqDistToLine2d(a, b);
}

template <class Line>
double distToLine2d ( const Line& a, const Line& b )
{
    return sqrt(sqDistToLine(a, b));
}
double distToLine ( const Line2d& a, const Line2d& b )
{
    return sqrt(distToLine2d(a, b));
}

// Measure the maximum distance between line segments.
// -- Arguments --
// a : line segment to compare.
// b : line segment to compare.
// -- Returns --
// The maximum square distance between the line segments. If one
// lies entirely on the other, returns 0.0.
template <class Line>
double sqMaxDistToLine2d ( const Line& a, const Line& b )
{
    double d = 0.0;
    
    // if zero-length segment, treat as a point
    if ( a.len < 0.000001 ) return sqDistToLine2d(a.a, b);
    
    // get distances between end-points
    double aaba = sqDist2d(a.a, b.a);
    double aabb = sqDist2d(a.a, b.b);
    double abba = sqDist2d(a.b, b.a);
    double abbb = sqDist2d(a.b, b.b);
    
    // check if line segments are parallel
    int isNotParallel = isRight(a.r, b.r);
    if ( !isNotParallel ) {
        // get distances between end-points and line segments
        double aab = sqDistToLine2d(a.a, b);
        double abb = sqDistToLine2d(a.b, b);
        double baa = sqDistToLine2d(b.a, a);
        double bba = sqDistToLine2d(b.b, a);
        
        // check if one line lies within the other
        if ( (aab < 0.0000001 && abb < 0.0000001)
                || (baa < 0.0000001 && bba < 0.0000001) )
            return 0.0;
    }
    
    // get largest sq distances between segment ends
    d = aaba;
    if ( aabb > d ) d = aabb;
    if ( abba > d ) d = abba;
    if ( abbb > d ) d = abbb;
    
    return d;
}
double sqMaxDistToLine ( const Line2d& a, const Line2d& b )
{
    return sqMaxDistToLine2d(a, b);
}

template <class Line>
double maxDistToLine2d ( const Line& a, const Line& b )
{
    return sqrt(sqMaxDistToLine(a, b));
}
double maxDistToLine ( const Line2d& a, const Line2d& b )
{
    return sqrt(sqMaxDistToLine2d(a, b));
}

// Measure the angle between two lines.
// -- Arguments --
// a : line to compare.
// b : line to compare.
// -- Returns --
// The angle between two lines in radians. Returns a negative value
// if an angle cannot be determined (i.e., if a or b is 0-length).
double angle2d ( const Line2d& a, const Line2d& b )
{
    double denom = a.len * b.len;
    if ( denom == 0.0 ) return -1;
    
    double num = a.dx * b.dx + a.dy * b.dy;
    return acos(num / denom);
};


#endif // YOUNG_GEOMETRY_LINE_20201211