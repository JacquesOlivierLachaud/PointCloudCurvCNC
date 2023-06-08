/**
 Copyright (c) 2021
 Jacques-Olivier Lachaud (\c jacques-olivier.lachaud@univ-savoie.fr)
 Laboratory of Mathematics (CNRS, UMR 5127), University Savoie Mont Blanc, France,

All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
* Redistributions of source code must retain the above copyright
  notice, this list of conditions and the following disclaimer.

* Redistributions in binary form must reproduce the above copyright
  notice, this list of conditions and the following disclaimer in the
  documentation and/or other materials provided with the distribution.

* Neither the name of the <organization> nor the names of its
  contributors may be used to endorse or promote products derived from
  this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT
HOLDER> BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

/**
 * @file RegularPolyhedron.h
 * @author Jacques-Olivier Lachaud (\c jacques-olivier.lachaud@univ-savoie.fr)
 * Laboratory of Mathematics (CNRS, UMR 5127), University Savoie Mont Blanc, France
 *
 * @date 2021/04/03
 *
 */

#pragma once

#include "Common.h"

/// Represents a regular 3D polyhedron with either 6 or 12 sides.
///
/// @tparam TRealPoint any model of 3D point
template <typename TRealPoint>
struct RegularPolyhedron
{
  typedef int        Index;
  typedef TRealPoint RealPoint;
  typedef RealPoint  RealVector;
  typedef double     Scalar;
  typedef std::vector< int > Face;
  
  Scalar                   _size; ///< the radius of the enclosing sphere
  Scalar                   _face_radius; ///< the radius of every face
  std::vector< RealPoint > _points; ///< the vertices of the regular polyhedron
  std::vector< Face >      _faces;  ///< the vertex indices for each face
  std::vector< Scalar >    _intercept; ///< the intercept of the plane supporting each face
  std::vector< RealPoint > _centers; ///< the centers of each face
  std::vector< RealPoint > _normals; ///< the normal of the plane supporting each face
  std::vector< RealPoint > _bases0;  ///< a tangent vector of the plane supporting each face
  std::vector< RealPoint > _bases1;  ///< another  tangent vector of the plane supporting each face

  /// Constructor
  /// @param N the number of sides, either 6 (cube) or 12 (dodecahedron).
  /// @param s the radius of the bounding sphere.
  RegularPolyhedron( int N = 6, Scalar s = 1.0 )
  {
    _size = s;
    switch (N) {
    case 6: initCube(); break;
    case 12: initDodecahedron(); break;
    }
    bool ok = isInside( RealPoint( 0.0, 0.0, 0.0 ) );
    if ( !ok )
      std::cerr << "Bad initialization: center is not inside." << std::endl;
  }

  /// @param[in] p any 3D point
  /// @return 'true' if and only if \a p is inside the polyhedron.
  bool isInside( RealPoint p ) const
  {
    for ( Index i = 0; i < _faces.size(); i++ )
      if ( _normals[ i ].dot( p ) - _intercept[ i ] > 0.000001 ) return false;
    return true;
  }

  /// @return a couple `(p,n)` where `p` is a point randomly chosen on
  /// the polyhedron surface and `n` is its normal vector.
  std::pair< RealPoint, RealVector >
  randomPointNormal() const
  {
    Index f = rand() % _faces.size();
    RealPoint p;
    do {
      double x = (2.0 * rand01() - 1.0) * _face_radius;
      double y = (2.0 * rand01() - 1.0) * _face_radius;
      p = _centers[ f ];
      p += _bases0[ f ] * x;
      p += _bases1[ f ] * y;
    } while ( ! isInside( p ) );
    return std::make_pair( p, _normals[ f ] );
  }

  // --------------------------------------------------------------------------

  /// Computes the normal vectors to each face.
  void computeNormals()
  {
    const Index nbf = _faces.size();
    _normals.resize( nbf );
    _intercept.resize( nbf );
    for ( Index f = 0; f < nbf; f++ )
      {
        RealVector n =
          ( _points[ _faces[ f ][ 1 ] ] - _points[ _faces[ f ][ 0 ] ] )
          .cross( _points[ _faces[ f ][ 2 ] ] - _points[ _faces[ f ][ 0 ] ] );
        n /= n.norm();
        Scalar    mu = n.dot( _points[ _faces[ f ][ 0 ] ] );
        bool      in = n.dot( RealPoint( 0.0, 0.0, 0.0 ) ) <= mu;
        _normals  [ f ] = in ?  n :  -n;
        _intercept[ f ] = in ? mu : -mu;
      }
  }

  /// Computes tangent vector bases for each face.  
  void computeBases()
  {
    const Index nbf = _faces.size();
    _bases0.resize( nbf );
    _bases1.resize( nbf );
    for ( Index f = 0; f < nbf; f++ )
      {
        auto u = _points[ _faces[ f ][ 1 ] ] - _points[ _faces[ f ][ 0 ] ];
        auto v = _normals[ f ].cross( u );
        _bases0[ f ] = u / u.norm();
        _bases1[ f ] = v / v.norm();
      }
  }

  /// Computes the centers of each face.
  void computeCenters()
  {
    const Index nbf = _faces.size();
    _centers.resize( nbf );
    for ( Index f = 0; f < nbf; f++ )
      {
        RealPoint c( 0.0, 0.0, 0.0 );
        for ( auto i : _faces[ f ] )
          c += _points[ i ];
        c /= _faces[ f ].size();
        _centers[ f ] = c;
      }
  }

  /// Initializes a cube as a regular polyhedron.
  void initCube()
  {
    Scalar s = _size;
    _points.push_back( RealPoint(  s,  s,  s ) );
    _points.push_back( RealPoint( -s,  s,  s ) );    
    _points.push_back( RealPoint(  s, -s,  s ) );
    _points.push_back( RealPoint( -s, -s,  s ) );    
    _points.push_back( RealPoint(  s,  s, -s ) );
    _points.push_back( RealPoint( -s,  s, -s ) );    
    _points.push_back( RealPoint(  s, -s, -s ) );
    _points.push_back( RealPoint( -s, -s, -s ) );
    _faces.push_back( Face { 0, 1, 3, 2 } );
    _faces.push_back( Face { 0, 2, 6, 4 } );
    _faces.push_back( Face { 0, 4, 5, 1 } );
    _faces.push_back( Face { 7, 5, 4, 6 } );
    _faces.push_back( Face { 7, 6, 2, 3 } );
    _faces.push_back( Face { 7, 3, 1, 5 } );
    computeNormals();
    computeBases();
    computeCenters();
    _face_radius = _size * sqrt(2);
  }

  /// Initializes a dodecahedron as a regular polyhedron.
  void initDodecahedron()
  {
    Scalar s = _size;
    Scalar h = 2.0 / (1.0+sqrt(5));
    _points.push_back( RealPoint(  s,  s,  s ) );
    _points.push_back( RealPoint( -s,  s,  s ) );    
    _points.push_back( RealPoint(  s, -s,  s ) );
    _points.push_back( RealPoint( -s, -s,  s ) );    
    _points.push_back( RealPoint(  s,  s, -s ) );
    _points.push_back( RealPoint( -s,  s, -s ) );    
    _points.push_back( RealPoint(  s, -s, -s ) );
    _points.push_back( RealPoint( -s, -s, -s ) );
    
    _points.push_back( RealPoint(  s*(1.0+h),  s*(1.0-h*h),  0.0 ) );
    _points.push_back( RealPoint(  s*(1.0+h), -s*(1.0-h*h),  0.0 ) );
    _points.push_back( RealPoint( -s*(1.0+h),  s*(1.0-h*h),  0.0 ) );
    _points.push_back( RealPoint( -s*(1.0+h), -s*(1.0-h*h),  0.0 ) );
    _points.push_back( RealPoint(  0.0,  s*(1.0+h),  s*(1.0-h*h) ) );
    _points.push_back( RealPoint(  0.0,  s*(1.0+h), -s*(1.0-h*h) ) );
    _points.push_back( RealPoint(  0.0, -s*(1.0+h),  s*(1.0-h*h) ) );
    _points.push_back( RealPoint(  0.0, -s*(1.0+h), -s*(1.0-h*h) ) );
    _points.push_back( RealPoint(  s*(1.0-h*h),  0.0,  s*(1.0+h) ) );
    _points.push_back( RealPoint( -s*(1.0-h*h),  0.0,  s*(1.0+h) ) );
    _points.push_back( RealPoint(  s*(1.0-h*h),  0.0, -s*(1.0+h) ) );
    _points.push_back( RealPoint( -s*(1.0-h*h),  0.0, -s*(1.0+h) ) );

    // faces associées à +x
    _faces.push_back( Face { 8, 9, 2, 16, 0 } );
    _faces.push_back( Face { 9, 8, 4, 18, 6 } );
    // faces associées à -x
    _faces.push_back( Face { 11, 10, 1, 17, 3 } );
    _faces.push_back( Face { 10, 11, 7, 19, 5 } );
    // faces associées à +y
    _faces.push_back( Face { 12, 13, 4, 8, 0 } );
    _faces.push_back( Face { 13, 12, 1, 10, 5 } );
    // faces associées à -y
    _faces.push_back( Face { 15, 14, 2, 9, 6 } );
    _faces.push_back( Face { 14, 15, 7, 11, 3 } );
    // faces associées à +z
    _faces.push_back( Face { 16, 17, 1, 12, 0 } );
    _faces.push_back( Face { 17, 16, 2, 14, 3 } );
    // faces associées à -z
    _faces.push_back( Face { 19, 18, 4, 13, 5 } );
    _faces.push_back( Face { 18, 19, 7, 15, 6 } );
    
    computeNormals();
    computeBases();
    computeCenters();
    _face_radius = _size;
  }

  /// @return a random floating point number between 0 and 1.
  static Scalar rand01()
  {
    return (double) rand() / (double) RAND_MAX;
  }
  
};
