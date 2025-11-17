/*******************************************************************************
 * TIN-based Real-time Integrated Basin Simulator (tRIBS)
 * Distributed Hydrologic Model
 *
 * Copyright (c) 2025. tRIBS Developers
 *
 * See LICENSE file in the project root for full license information.
 ******************************************************************************/

/***************************************************************************
**
**  geometry.h:  Definitions of some simple geometry classes.
**
***************************************************************************/

#ifndef GEOMETRY_H
#define GEOMETRY_H

//=========================================================================
//
//
//                  Section 1: geometry Class Declarations
//
//
//=========================================================================

class Point2D
{
  public:
   Point2D();
   Point2D( double, double );
   const Point2D &operator=( const Point2D & );
   double x;
   double y;
};

class Point3D
{
  public:
   Point3D();
   Point3D( double, double, double );
   const Point3D &operator=( const Point3D & );
   double x;
   double y;
   double z;
};

//=========================================================================
//
//
//                  Section 2: geometry Inline Functions
//
//
//=========================================================================

inline Point2D::Point2D() {
   x = y = 0.0;
}

inline Point2D::Point2D( double ix, double iy ){
   x = ix;
   y = iy;
}

inline const Point2D &Point2D::operator=(const Point2D &right ){
    x = right.x;
    y = right.y;
    return *this;
}

inline Point3D::Point3D() {
   x = y = z = 0.0;
}

inline Point3D::Point3D( double ix, double iy, double iz ){
   x = ix;
   y = iy;
   z = iz;
}

inline const Point3D &Point3D::operator=(const Point3D &right ){
    x = right.x;
    y = right.y;
    z = right.z;
    return *this;
}

#endif

//=========================================================================
//
//
//                         End of geometry.h
//
//
//=========================================================================
