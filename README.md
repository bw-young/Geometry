# Geometry
Objects and operations for representing geometric objects and problems.

This is an ongoing work in progress and is liable to change at any time without warning. For example, a planned changes includes generalization of the objects to n dimensions, which will invalidate the present Point2d, Point3d, etc., etc. framework currently in place.

Likewise, this documentation will expand and become more complex as needed.

## Objects

### Point2d
Header file: Point.hpp

```cpp
class Point2d {double x, y;};

Point2d p (x,y);  // construct with (x,y); default (0,0)
Point2d p2 (p);   // copy constructor

p2 = p;           // assignment copies the object's contents
p < p2;           // less than if x is smaller or, if x is the same, if y is smaller.
p == p2;          // check if x and y are exactly the same.
p != p2;          // check if x and y are exactly the same.
```

### Point3d
Header file: Point.hpp

```cpp
class Point3d {double x, y, z;};

Point2d p (x,y,z);  // construct with (x,y,z); default (0,0,0)
Point2d p2 (p);     // copy constructor

p2 = p;             // assignment copies the object's contents
p < p2;             // true if x in p is smaller than x in p2 or, if x is the same, if y is smaller.
p == p2;            // true if x and y are exactly the same.
p != p2;            // !(p == p2)
```

### Line
Header file: Line.hpp

```cpp
class Line2d{
  Point2d a;      // segment start point
  Point2d b;      // segment end point
  Point2d r;      // segment direction vector a->b
  double dx;      // b.x - a.x
  double dy;      // b.y - a.y
  double len;     // segment length
}; // Line2d

Line2d s (p1,p2); // construct with (p1,p2) and computes all other members of s; default (Point2d(),Point2d())
Line2d s2 (s);    // copy constructor
s2 = s;           // assignment copies the object's contents
s < s2;           // true if both points in s are less than either point in s2
s == s2;          // true if points in s are equal to the same points in s2
s != s2;          // !(s == s2)
```

## Operations
```cpp
// square euclidean distance between points a and b
template <class Point> double sqDist2d (Point a, Point b);
double sqDist(Point2d a, Point2d b);
  
// euclidean distance between points a and b
template <class Point> double dist2d (Point a, Point b);
double dist(Point2d a, Point2d b);

// euclidean distance to ray r0->r1 intersection with line segment p0->p1, negative otherwise
template <class Point> double raySegment_intersect2d (Point r0, Point r1, Point p0, Point p1);

// square euclidean distance between point a and line cd
template <class Point, class Line> double sqDistToLine2d (Point a, Line cd);
double sqDistToLine (Point2d a, Line2d cd);

// euclidean distance between point a and line cd
template <class Point, class Line> double distToLine2d (Point a, Line cd);
double distToLine (Point2d a, Line2d cd);

// square euclidean distance between lines ab and cd
template <class Line> double sqDistToLine2d (Line ab, Line cd);
double sqDistToLine (Line2d ab, Line2d cd);

// euclidean distance between lines ab and cd
template <class Line> double distToLine2d (Line ab, Line cd);
double distToLine (Line2d ab, Line2d cd);
```
