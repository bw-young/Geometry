/////////////////////////////////////////////////////////////////////
// Objects and methods for managing geometric points.              //
/////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////
// -- HISTORY ---------------------------------------------------- //
// 04/23/2020 - Brennan Young                                      //
// - created                                                       //
// - migrated Point2d, Point3d, inPoly, raySegment_intersect2d,    //
//   and rayPoly_intersect2d from GridObjPolyOps.h.                //
// 04/27/2020 - Brennan Young                                      //
// - Point2d and Point3d assignment operator now returns *this.    //
// 12/11/2020                                                      //
// - migrated Point2d, Point3d from ObjectPolyOps.h.               //
// - added sqDist2d, sqDist, dist2d, and dist.                     //
/////////////////////////////////////////////////////////////////////

#ifndef YOUNG_GEOMETRY_POINT_20201211
#define YOUNG_GEOMETRY_POINT_20201211


// Basic 2d point object.
class Point2d {
public:
    double x, y;
    
    // constructors, destructor
    Point2d ( double X=0.0, double Y=0.0 ) : x(X), y(Y) {}
    Point2d ( const Point2d& p ) : x(p.x), y(p.y) {}
    ~Point2d () {}
    
    // operators
    Point2d& operator= ( const Point2d& p )
    {
        if ( &p == this ) return *this;
        x = p.x;
        y = p.y;
        return *this;
    }
    bool operator< ( const Point2d& p ) const
    {
        // sort by y, then x
        if ( y < p.y ) return true;
        if ( y > p.y ) return false;
        return x < p.x;
    }
    bool operator== ( const Point2d& p ) const
    {
        return x == p.x && y == p.y;
    }
    bool operator!= ( const Point2d& p ) const
    {
        return !(x == p.x && y == p.y);
    }
}; // Point2d


// Basic 3D point object.
class Point3d {
public:
    double x, y, z;
    
    // constructors, destructor
    Point3d ( double X=0.0, double Y=0.0, double Z=0.0 ) :
        x(X), y(Y),z(Z) {}
    Point3d ( const Point3d& p ) : x(p.x), y(p.y), z(p.z) {}
    ~Point3d () {}
    
    // operators
    Point3d& operator= ( const Point3d& p )
    {
        if ( &p == this ) return *this;
        x = p.x;
        y = p.y;
        z = p.z;
        return *this;
    }
    bool operator< ( const Point3d& p ) const
    {
        // sort by z, then y, then x
        if ( z < p.z ) return true;
        if ( z > p.z ) return false;
        if ( y < p.y ) return true;
        if ( y > p.y ) return false;

        return x < p.x;
    }
    bool operator== ( const Point3d& p ) const
    {
        return x == p.x && y == p.y && z == p.z;
    }
    bool operator!= ( const Point3d& p ) const
    {
        return !(x == p.x && y == p.y && z != p.z);
    }
}; // Point3d


/////////////////////////////////////////////////////////////////////
// Primitive Operations                                            //
/////////////////////////////////////////////////////////////////////


// Get the distance between two point objects.
template <class Point>
double sqDist2d ( const Point& a, const Point& b )
{
    double dx = b.x - a.x;
    double dy = b.y - a.y;
    return dx*dx + dy*dy;
}
double sqDist ( const Point2d& a, const Point2d& b )
{
    return sqDist2d(a, b);
}

template <class Point>
double dist2d ( const Point& a, const Point& b )
{
    return sqrt(sqDist2d(a, b));
}
double dist ( const Point2d& a, const Point2d& b )
{
    return dist2d(a, b);
}


#endif // YOUNG_GEOMETRY_POINT_20201211