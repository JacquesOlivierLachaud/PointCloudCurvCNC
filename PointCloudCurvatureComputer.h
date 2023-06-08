/**
 Copyright (c) 2022
 Jacques-Olivier Lachaud (\c jacques-olivier.lachaud@univ-savoie.fr)
 Laboratory of Mathematics (CNRS, UMR 5807), University of Savoie, France,
 
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

#pragma once

#include "Common.h"
#include "LinearKDTree.h"
#include "CorrectedNormalCurrentFormulaEigen.h"

///////////////////////////////////////////////////////////////////////////////
// HEXAGRAM TRIANGLE GENERATION
///////////////////////////////////////////////////////////////////////////////

/// Structure to generate the (two) local triangles in case of Hexagram method.
struct HexagramGeneration {
  typedef Eigen::Vector3d    RealPoint;
  typedef Eigen::Vector3d    RealVector;
  typedef std::size_t        Index;
  typedef std::vector<Index> Indices;
  typedef double             Scalar;

  Scalar    _avgnormals;
  std::array< Index,   6 > _vertices;
  std::array< Scalar,    6 > _distance2;
  std::array< RealPoint, 6 > _targets;
  std::array< Scalar,    6 > _cos;
  std::array< Scalar,    6 > _sin;
  
  HexagramGeneration( Scalar average_normals )
    : _avgnormals( average_normals )
  {
    for ( auto j = 0; j < 6; j++ )
      {
        const Scalar a = j * M_PI / 3.0;
        _cos[ j ] = cos( a );
        _sin[ j ] = sin( a );
      }
  }
  
  void bin( Index i,
            const Indices& indices,
            const std::vector<RealPoint>& points,
            const std::vector<RealPoint>& normals )
  {
    // Compute normal and maximum distance.
    const RealPoint c = points[ i ];
    Scalar  avgd = 0.0;
    RealVector n = normals[ i ];
    RealVector a( 0.0, 0.0, 0.0 );
    for ( auto j : indices )
      {
        avgd += ( points [ j ] - c ).norm();
        a    += normals[ j ];
      }
    a    /= a.norm();
    n     = ( 1.0 - _avgnormals ) * n + _avgnormals * a;
    n    /= n.norm();
    avgd /= indices.size();
    
    // Define basis for sector analysis.
    const int  m = ( fabs( n[ 0 ] ) > fabs( n[ 1 ] ) )
      ? ( ( fabs( n[ 0 ] ) > fabs( n[ 2 ] ) ) ? 0 : 2 )
      : ( ( fabs( n[ 1 ] ) > fabs( n[ 2 ] ) ) ? 1 : 2 );
    const RealVector e =
      ( m == 0 ) ? RealVector( 0.0, 1.0, 0.0  ) :
      ( m == 1 ) ? RealVector( 0.0, 0.0, 1.0  ) :
      RealVector( 1.0, 0.0, 0.0  );
    RealVector u = n.cross( e );
    RealVector v = n.cross( u );
    u /= u.norm();
    v /= v.norm();
    for ( auto k = 0; k < 6; k++ )
      {
        _vertices [ k ] = i;
        _distance2[ k ] = avgd * avgd;
        _targets   [ k ] = avgd * ( u * _cos[ k ] + v * _sin[ k ] );
      }
    // Compute closest points.
    for ( auto j : indices )
      {
        if ( i == j ) continue;
        const RealVector w = ( points[ j ] - c );
        for ( auto k = 0; k < 6; k++ )
          {
            const Scalar d2 = ( w - _targets[ k ] ).squaredNorm();
            if ( d2 < _distance2[ k ] ) {
              _vertices [ k ] = j;
              _distance2[ k ] = d2;
            }
          }
      }
  }

  /// Computes and returns 2 random triangles.
  void getTriangles( std::vector< std::size_t >& t_j,
                     std::vector< std::size_t >& t_k,
                     std::vector< std::size_t >& t_l )
  {
    t_j.push_back( _vertices[ 0 ] );
    t_j.push_back( _vertices[ 1 ] );    
    t_k.push_back( _vertices[ 2 ] );
    t_k.push_back( _vertices[ 3 ] );    
    t_l.push_back( _vertices[ 4 ] );
    t_l.push_back( _vertices[ 5 ] );    
  }

  /// Computes two triangles. Vectors are preallocated.
  /// @return the number of triangles put into the vectors, i.e. 2.
  std::size_t fillTriangles( std::vector< std::size_t >& t_j,
                             std::vector< std::size_t >& t_k,
                             std::vector< std::size_t >& t_l )
  {
    t_j[ 0 ] = _vertices[ 0 ];
    t_j[ 1 ] = _vertices[ 1 ];    
    t_k[ 0 ] = _vertices[ 2 ];
    t_k[ 1 ] = _vertices[ 3 ];    
    t_l[ 0 ] = _vertices[ 4 ];
    t_l[ 1 ] = _vertices[ 5 ];
    return 2;
  }

  /// All in one method that builds the two triangles
  std::size_t
  fillPointNormals( std::vector< RealVector >& pn,
                    Index i,
                    const Indices& indices,
                    const std::vector<RealPoint>& points,
                    const std::vector<RealPoint>& normals )
  {
    // Compute normal and maximum distance.
    const RealPoint c = points[ i ];
    Scalar  avgd = 0.0;
    const RealVector ni = normals[ i ];
    RealVector       n  = ni;
    RealVector a( 0.0, 0.0, 0.0 );
    for ( auto j : indices )
      {
        avgd += ( points [ j ] - c ).norm();
        a    += normals[ j ];
      }
    a    /= a.norm();
    n     = ( 1.0 - _avgnormals ) * n + _avgnormals * a;
    n    /= n.norm();
    avgd /= indices.size();
    
    // Define basis for sector analysis.
    const int  m = ( fabs( n[ 0 ] ) > fabs( n[ 1 ] ) )
      ? ( ( fabs( n[ 0 ] ) > fabs( n[ 2 ] ) ) ? 0 : 2 )
      : ( ( fabs( n[ 1 ] ) > fabs( n[ 2 ] ) ) ? 1 : 2 );
    const RealVector e =
      ( m == 0 ) ? RealVector( 0.0, 1.0, 0.0  ) :
      ( m == 1 ) ? RealVector( 0.0, 0.0, 1.0  ) :
      RealVector( 1.0, 0.0, 0.0  );
    RealVector u = n.cross( e );
    RealVector v = n.cross( u );
    u /= u.norm();
    v /= v.norm();
    for ( auto k = 0; k < 6; k++ )
      {
        _vertices [ k ] = i;
        _distance2[ k ] = avgd * avgd;
        _targets  [ k ] = avgd * ( u * _cos[ k ] + v * _sin[ k ] );
      }
    const RealVector gn = _targets[ 0 ].cross( _targets[ 1 ] );
    if ( n.dot( gn ) < 0.0 ) n *= -1.0;
    // Compute closest points.
    for ( auto j : indices )
      {
        if ( i == j ) continue;
        const RealVector w = ( points[ j ] - c );
        for ( auto k = 0; k < 6; k++ )
          {
            const Scalar d2 = ( w - _targets[ k ] ).squaredNorm();
            if ( d2 < _distance2[ k ] ) {
              _vertices [ k ] = j;
              _distance2[ k ] = d2;
            }
          }
      }
    auto data = pn.data();
    const auto j1 = _vertices[ 0 ];
    const auto j2 = _vertices[ 1 ];    
    const auto k1 = _vertices[ 2 ];
    const auto k2 = _vertices[ 3 ];    
    const auto l1 = _vertices[ 4 ];
    const auto l2 = _vertices[ 5 ];
    *data++ = points[ j1 ];
    *data++ = points[ k1 ];
    *data++ = points[ l1 ];
    *data++ = normals[ j1 ];
    *data++ = normals[ k1 ];
    *data++ = normals[ l1 ];            
    *data++ = points[ j2 ];
    *data++ = points[ k2 ];
    *data++ = points[ l2 ];
    *data++ = normals[ j2 ];
    *data++ = normals[ k2 ];
    *data++ = normals[ l2 ];
    return 2;
  }
  
};


///////////////////////////////////////////////////////////////////////////////
// AVERAGED HEXAGRAM GENERATION
///////////////////////////////////////////////////////////////////////////////

/// Structure to generate the (two) local triangles in case of Averaged Hexagram method.
struct AvgHexagramGeneration {
  typedef Eigen::Vector3d    RealPoint;
  typedef Eigen::Vector3d    RealVector;
  typedef std::size_t        Index;
  typedef std::vector<Index> Indices;
  typedef double             Scalar;

  Scalar    _avg_normals_coef;
  std::array< RealPoint, 6 > _targets;
  std::array< RealVector,6 > _avg_normals;
  std::array< RealVector,6 > _avg_points;  
  std::array< Scalar,    6 > _cos;
  std::array< Scalar,    6 > _sin;
  std::array< size_t,    6 > _nb;
  
  AvgHexagramGeneration( Scalar average_normals )
    : _avg_normals_coef( average_normals )
  {
    for ( auto j = 0; j < 6; j++ )
      {
        const Scalar a = j * M_PI / 3.0;
        _cos[ j ] = cos( a );
        _sin[ j ] = sin( a );
      }
  }
  
  void bin( Index i,
            const Indices& indices,
            const std::vector<RealPoint>& points,
            const std::vector<RealPoint>& normals )
  {
    // Compute normal and maximum distance.
    const RealPoint c = points[ i ];
    Scalar  avgd = 0.0;
    RealVector n = normals[ i ];
    RealVector a( 0.0, 0.0, 0.0 );
    for ( auto j : indices )
      {
        avgd += ( points [ j ] - c ).norm();
        a    += normals[ j ];
      }
    a    /= a.norm();
    n     = ( 1.0 - _avg_normals_coef ) * n + _avg_normals_coef * a;
    n    /= n.norm();
    avgd /= indices.size();
    
    // Define basis for sector analysis.
    const int  m = ( fabs( n[ 0 ] ) > fabs( n[ 1 ] ) )
      ? ( ( fabs( n[ 0 ] ) > fabs( n[ 2 ] ) ) ? 0 : 2 )
      : ( ( fabs( n[ 1 ] ) > fabs( n[ 2 ] ) ) ? 1 : 2 );
    const RealVector e =
      ( m == 0 ) ? RealVector( 0.0, 1.0, 0.0  ) :
      ( m == 1 ) ? RealVector( 0.0, 0.0, 1.0  ) :
      RealVector( 1.0, 0.0, 0.0  );
    RealVector u = n.cross( e );
    RealVector v = n.cross( u );
    u /= u.norm();
    v /= v.norm();
    const RealVector zero( 0.0, 0.0, 0.0 );
    for ( auto k = 0; k < 6; k++ )
      {
        _targets    [ k ] = avgd * ( u * _cos[ k ] + v * _sin[ k ] );
        _avg_normals[ k ] = zero;
        _avg_points [ k ] = zero;
        _nb         [ k ] = 0;
      }
    // Compute closest points.
    for ( auto n = 0; n < indices.size(); n++ )
      {
        const auto       j = indices[ n ];
        const RealVector w = ( points[ j ] - c );
        auto best_k  = 0;
        auto best_d2 = ( w - _targets[ 0 ] ).squaredNorm();
        for ( auto k = 1; k < 6; k++ )
          {
            const Scalar d2 = ( w - _targets[ k ] ).squaredNorm();
            if ( d2 < best_d2 ) {
              best_k  = k;
              best_d2 = d2;
            }
          }
        _avg_normals[ best_k ] += normals[ j ];
        _avg_points [ best_k ] += points [ j ];
        _nb         [ best_k ] += 1;
      }
    for ( auto k = 0; k < 6; k++ )
      {
        if ( _nb[ k ] == 0 )
          {
            _avg_normals[ k ] = normals[ i ];
            _avg_points [ k ] = points [ i ];
          }
        else
          {
            _avg_normals[ k ] /= _avg_normals[ k ].norm();
            _avg_points [ k ] /= _nb[ k ];
          }
      }
  }

  /// All in one method that builds the two triangles.
  std::size_t
  fillPointNormals( std::vector< RealVector >& pn,
                    Index i,
                    const Indices& indices,
                    const std::vector<RealPoint>& points,
                    const std::vector<RealPoint>& normals )
  {
    // Compute normal and maximum distance.
    const RealPoint c = points[ i ];
    Scalar  avgd = 0.0;
    const RealVector ni = normals[ i ];
    RealVector       n  = ni;
    RealVector a( 0.0, 0.0, 0.0 );
    for ( auto j : indices )
      {
        avgd += ( points [ j ] - c ).norm();
        a    += normals[ j ];
      }
    a    /= a.norm();
    n     = ( 1.0 - _avg_normals_coef ) * n + _avg_normals_coef * a;
    n    /= n.norm();
    avgd /= indices.size();
    
    // std::cout << "Normal is " << n << std::endl;
    // Define basis for sector analysis.
    const int  m = ( fabs( n[ 0 ] ) > fabs( n[ 1 ] ) )
      ? ( ( fabs( n[ 0 ] ) > fabs( n[ 2 ] ) ) ? 0 : 2 )
      : ( ( fabs( n[ 1 ] ) > fabs( n[ 2 ] ) ) ? 1 : 2 );
    const RealVector e =
      ( m == 0 ) ? RealVector( 0.0, 1.0, 0.0  ) :
      ( m == 1 ) ? RealVector( 0.0, 0.0, 1.0  ) :
      RealVector( 1.0, 0.0, 0.0  );
    RealVector u = n.cross( e );
    RealVector v = n.cross( u );
    u /= u.norm();
    v /= v.norm();
    const RealVector zero( 0.0, 0.0, 0.0 );
    for ( auto k = 0; k < 6; k++ )
      {
        _targets    [ k ] = avgd * ( u * _cos[ k ] + v * _sin[ k ] );
        _avg_normals[ k ] = zero;
        _avg_points [ k ] = zero;
        _nb         [ k ] = 0;
      }
    // Compute closest points.
    for ( auto m = 0; m < indices.size(); m++ )
      {
        const auto       j = indices[ m ];
        const RealVector w = ( points[ j ] - c );
        auto best_k  = 0;
        auto best_d2 = ( w - _targets[ 0 ] ).squaredNorm();
        for ( auto k = 1; k < 6; k++ )
          {
            const Scalar d2 = ( w - _targets[ k ] ).squaredNorm();
            if ( d2 < best_d2 ) {
              best_k  = k;
              best_d2 = d2;
            }
          }
        _avg_normals[ best_k ] += normals[ j ];
        _avg_points [ best_k ] += points [ j ];
        _nb         [ best_k ] += 1;
      }
    auto data = pn.data();
    *data++   = ( _nb[ 0 ] == 0 ) ? points [ i ] : _avg_points [ 0 ] / _nb[ 0 ];
    *data++   = ( _nb[ 2 ] == 0 ) ? points [ i ] : _avg_points [ 2 ] / _nb[ 2 ];
    *data++   = ( _nb[ 4 ] == 0 ) ? points [ i ] : _avg_points [ 4 ] / _nb[ 4 ];
    *data++   = ( _nb[ 0 ] == 0 ) ? normals[ i ] : _avg_normals[ 0 ] / _avg_normals[ 0 ].norm();
    *data++   = ( _nb[ 2 ] == 0 ) ? normals[ i ] : _avg_normals[ 2 ] / _avg_normals[ 2 ].norm();
    *data++   = ( _nb[ 4 ] == 0 ) ? normals[ i ] : _avg_normals[ 4 ] / _avg_normals[ 4 ].norm();
    *data++   = ( _nb[ 1 ] == 0 ) ? points [ i ] : _avg_points [ 1 ] / _nb[ 1 ];
    *data++   = ( _nb[ 3 ] == 0 ) ? points [ i ] : _avg_points [ 3 ] / _nb[ 3 ];
    *data++   = ( _nb[ 5 ] == 0 ) ? points [ i ] : _avg_points [ 5 ] / _nb[ 5 ];
    *data++   = ( _nb[ 1 ] == 0 ) ? normals[ i ] : _avg_normals[ 1 ] / _avg_normals[ 1 ].norm();
    *data++   = ( _nb[ 3 ] == 0 ) ? normals[ i ] : _avg_normals[ 3 ] / _avg_normals[ 3 ].norm();
    *data++   = ( _nb[ 5 ] == 0 ) ? normals[ i ] : _avg_normals[ 5 ] / _avg_normals[ 5 ].norm();
    return 2;
  }
   
};


///////////////////////////////////////////////////////////////////////////////
// Main class for computing curvature measures on point clouds.
///////////////////////////////////////////////////////////////////////////////

/// Main class for computing curvature measures on point clouds.
struct PointCloudCurvatureComputer {
  typedef Eigen::VectorXd  DenseVector;
  typedef Eigen::MatrixXd  DenseMatrix;
  typedef Eigen::Vector3d  RealPoint;
  typedef Eigen::Vector3d  RealVector;
  typedef Eigen::Matrix3d  RealTensor;
  typedef LinearKDTree< RealPoint, 3 > PointKDTree;
  typedef typename PointKDTree::Indices Indices;
  typedef std::size_t Size;
  typedef double      Scalar;

  enum class Method {
    UniformGeneration,
    IndependentGeneration,
    HexagramGeneration,
    AvgHexagramGeneration    
  };
  
  std::random_device _rd; ///< the random device
  std::mt19937       _rg; ///< the pseudo-random number generator.
  const PointKDTree* _ptr_tree; ///< the k-d-tree containing all the points.
  Size               _knn;      ///< the targetted number of neighbors
  Scalar             _best_rho; ///< ball radius that contains approx _k neighbors.
  Method             _method;   ///< methods for generating random triangles
  Size               _maxtriangles; ///< maximum number of triangles used for computation
  Scalar             _avgnormals; ///< weight of average normals around for the binner.

  /// Default constructor.
  /// The object is invalid and should be initialized with \ref init.
  PointCloudCurvatureComputer()
    : _rd(), _rg( _rd() ), _ptr_tree( nullptr ), _knn( 0 ), _best_rho( 0.0 ) {}

  /// Initializes the object with the range of points organized in a
  /// k-d-tree structure \a tree as well as the targeted number of
  /// neighbors \a  for each curvature computation.
  ///
  /// @param tree the k-d-tree structure containing all the points.
  /// @param knn the targeted number of neighbors for each point,
  /// which are used for curvature computations.
  void init( const PointKDTree& tree, int knn = 25 )
  {
    _ptr_tree = &tree;
    _knn      = knn;
    _best_rho = _ptr_tree->findRadius( knn, 0 );
    resetParameters();
  }

  /// @return the current number of points.
  Size nbPoint() const
  {
    return ( _ptr_tree == nullptr ) ? 0 : _ptr_tree->size();
  }

  /// @return the estimated radius such that every point has
  /// approximately _knn neighbors within the ball with such radius.
  Scalar bestRadius() const
  {
    return _best_rho;
  }
  
  void resetParameters()
  {
    _method             = Method::UniformGeneration;
    _maxtriangles       = 100; ///< maximum number of triangles used for computation
  }

  
  /// Choose the sampling method that builds triangles in the most
  /// uniformly random way, just picking triangle vertices in the
  /// neighborhood of the point of evaluattion in an arbitrary way.
  ///
  /// @param maxtriangles the maximum number of triangles used for computation
  void chooseUniformGenerationMethod( int  maxtriangles )
  {
    _method             = Method::UniformGeneration;
    _maxtriangles       = maxtriangles;
  }

  /// Choose the sampling method that tries to build independent
  /// random triangles without common vertices.
  ///
  /// @param maxtriangles the maximum number of generated triangles
  /// (here there at most `_knn / 3` generated triangles.
  void chooseIndependentGenerationMethod( int  maxtriangles )
  {
    _method             = Method::IndependentGeneration;
    _maxtriangles       = maxtriangles;
  }

  /// Choose the method that builds two triangles that surround the
  /// point of evaluation, approximately equilateral.
  ///
  /// @param average_normals_weight the weight used to balance between
  /// the normal to the current point of interest and the normals of
  /// its neighbors.
  void chooseHexagramGenerationMethod( Scalar average_normals_weight = 0.5 )
  {
    _method             = Method::HexagramGeneration;
    _avgnormals         = average_normals_weight;
  }

  /// Choose the method that builds two triangles that surround the
  /// point of evaluation, approximately equilateral, and choose the
  /// median normal estimation.
  ///
  /// @param average_normals_weight the weight used to balance between
  /// the normal to the current point of interest and the normals of
  /// its neighbors.
  void chooseAvgHexagramGenerationMethod( Scalar average_normals_weight = 0.5 )
  {
    _method             = Method::AvgHexagramGeneration;
    _avgnormals         = average_normals_weight;
  }
  
  
  /// Main method to compute the curvature measures of every point in the current
  /// k-d-tree of points. 
  ///
  /// @param[in] point_normals the vector giving for each point its normal vector.
  ///
  /// @return a tuple of vectors of size `size()`, in the following
  /// order: <1> area measure \f$ \mu^0 \f$; <2> mean curvature
  /// measure \f$ \mu^1 \f$; <3> Gauss curvature measure \f$ \mu^2
  /// \f$; <4,5,6,7,8,9> anisotropic curvature tensor measure \f$
  /// \mu^{XY} \f$, as the 6 vectors T11, T12, T13, T22, T23, T33;
  /// <10> timings.
  std::tuple< DenseVector, DenseVector, DenseVector,
              DenseVector, DenseVector, DenseVector, 
              DenseVector, DenseVector, DenseVector,
              std::vector<double> >
  computeCurvatureMeasures( const std::vector< RealVector >  & point_normals )
  {
    const auto nb = _ptr_tree->_points.size();
    _maxtriangles = std::max( _maxtriangles, 2UL );
    DenseVector A  ( nb ); // area
    DenseVector H  ( nb ); // mean curvatures
    DenseVector G  ( nb ); // Gaussian curvatures
    DenseVector T11( nb ); // curvatures tensor 11
    DenseVector T12( nb ); // curvatures tensor 12
    DenseVector T13( nb ); // curvatures tensor 13
    DenseVector T22( nb ); // curvatures tensor 22
    DenseVector T23( nb ); // curvatures tensor 23
    DenseVector T33( nb ); // curvatures tensor 33
    
    // Timing after the look-up
    std::vector<double> timings;
    // Store triangles
    std::vector< size_t >    t_j( _maxtriangles );
    std::vector< size_t >    t_k( _maxtriangles );
    std::vector< size_t >    t_l( _maxtriangles );    
    HexagramGeneration    STG ( _avgnormals );
    AvgHexagramGeneration SMTG( _avgnormals );
    bool smtg = ( _method == Method::AvgHexagramGeneration );
    // Compute curvature measures for each point.
    auto counter = 0;
    for ( auto i : _ptr_tree->_indices )
      {
	if ( counter >= 499499 ) std::cout << " " << counter << " " << i << std::endl;
	counter += 1;
	if ( i >= nb ) std::cerr << "Error" << std::endl;
        // Compute nearest neighbors and common stuff.
        const RealPoint  p = _ptr_tree->_points[ i ];
        auto nearest = _ptr_tree->kNeighborsAtLeast( p, _knn, _best_rho, true );
	if ( counter >= 499500 ) std::cout << "#nn=" << nearest.size() << std::endl;
        nearest.resize( _knn );
	if ( counter >= 499500 ) std::cout << "#nn=" << nearest.size() << std::endl;
        
        std::chrono::high_resolution_clock::time_point
	  t1 = std::chrono::high_resolution_clock::now();

        // Generate triangles
        auto nb_vt = 0; // number of valid generated triangles.
        if ( _method == Method::IndependentGeneration )
          {
            std::vector< size_t > indices( nearest.size() );
            for ( auto i = 0; i < indices.size(); i++ ) indices[ i ] = i;
            std::shuffle( indices.begin(), indices.end(), _rg );
            const auto maxt = std::min( _maxtriangles, indices.size() / 3 );
            auto i = 0;
            for ( ; nb_vt < maxt; ++nb_vt )
              {
                t_j[ nb_vt ] = nearest[ indices[ i++ ] ];
                t_k[ nb_vt ] = nearest[ indices[ i++ ] ];
                t_l[ nb_vt ] = nearest[ indices[ i++ ] ];                                
              }
          }
        else if ( _method == Method::HexagramGeneration )
          {
            STG.bin( i, nearest, _ptr_tree->_points, point_normals );
            nb_vt = STG.fillTriangles( t_j, t_k, t_l );              
          }
        else if ( _method == Method::AvgHexagramGeneration )
          {
            SMTG.bin( i, nearest, _ptr_tree->_points, point_normals );
            // nb_vt = SMTG.fillTriangles( t_j, t_k, t_l );
            nb_vt = 2;
          }
        if ( _method == Method::UniformGeneration
             || ( nb_vt == 0 ) )
          {
            for ( int t = 0; t < _maxtriangles; ++t )
              {
                auto j = nearest[ rand() % nearest.size() ];
                auto k = nearest[ rand() % nearest.size() ];
                auto l = nearest[ rand() % nearest.size() ];
                if ( j == k || j == l || k == l ) continue;
                t_j[ nb_vt ] = j;
                t_k[ nb_vt ] = k;
                t_l[ nb_vt ] = l;
                nb_vt       += 1;
              }
          }
	if ( counter >= 499500 ) std::cout << "nb_vt=" << nb_vt << std::endl;

        if ( nb_vt == 0 ) std::cerr << "No triangles." << std::endl;

        // Sum measures of all triangles and fix their orientation at the same time.
        double localA = 0.0;
        double localH = 0.0;
        double localG = 0.0;
        RealTensor localT = RealTensor::Zero();

        for ( auto t = 0; t < nb_vt; ++t )
          {
            const auto   j = t_j[ t ];
            const auto   k = t_k[ t ];
            const auto   l = t_l[ t ];

	    if ( counter >= 499500 )
	      std::cout << "j k l" << j << " " << k << " " << l << std::endl;
	    
            const auto& pj = smtg ? SMTG._avg_points [ t   ] : _ptr_tree->_points[ j ];
            const auto& pk = smtg ? SMTG._avg_points [ t+2 ] : _ptr_tree->_points[ k ];
            const auto& pl = smtg ? SMTG._avg_points [ t+4 ] : _ptr_tree->_points[ l ]; 
            const auto& nj = smtg ? SMTG._avg_normals[ t   ] : point_normals[ j ];
            const auto& nk = smtg ? SMTG._avg_normals[ t+2 ] : point_normals[ k ];
            const auto& nl = smtg ? SMTG._avg_normals[ t+4 ] : point_normals[ l ];
	    if ( counter >= 499500 ) {
	      std::cout << "pj pk pl " << pj << " " << pk << " " << pl << std::endl;
	      std::cout << "nj nk nl " << nj << " " << nk << " " << nl << std::endl;
	    }
        
            const auto  tA = CNCEigen::mu0InterpolatedU( pj, pk, pl, nj, nk, nl );
	    if ( counter >= 499500 ) std::cout << "tA = " << tA << std::endl; 
            if ( tA < -CNCEigen::epsilon )
              {
                localA -= tA;
                localH += CNCEigen::mu1InterpolatedU ( pj, pl, pk, nj, nl, nk );
                localG += CNCEigen::mu2InterpolatedU ( pj, pl, pk, nj, nl, nk );
                localT += CNCEigen::muXYInterpolatedU( pj, pl, pk, nj, nl, nk );
              }
            else if ( tA > CNCEigen::epsilon )
              { 
                localA += tA;
                localH += CNCEigen::mu1InterpolatedU ( pj, pk, pl, nj, nk, nl );
                localG += CNCEigen::mu2InterpolatedU ( pj, pk, pl, nj, nk, nl );
                localT += CNCEigen::muXYInterpolatedU( pj, pk, pl, nj, nk, nl );
              }
          } // for ( auto t = 0; t < t_j.size(); ++t )
        
        A  [ i ] = localA;
        H  [ i ] = localH;
        G  [ i ] = localG;
        T11[ i ] = localT( 0, 0 );
        T12[ i ] = 0.5 * ( localT( 0, 1 ) + localT( 1, 0 ) );
        T13[ i ] = 0.5 * ( localT( 0, 2 ) + localT( 2, 0 ) );
        T22[ i ] = localT( 1, 1 );
        T23[ i ] = 0.5 * ( localT( 1, 2 ) + localT( 2, 1 ) );
        T33[ i ] = localT( 2, 2 );
        
	if ( counter >= 499500 ) std::cout << "Filled vectors at " << i << std::endl;         
        std::chrono::high_resolution_clock::time_point
	  t2 = std::chrono::high_resolution_clock::now();
        timings.push_back
	  ( std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count() );
      } // for ( auto i : _ptr_tree->_indices )
    return std::make_tuple( A, H, G, T11, T12, T13, T22, T23, T33 , timings);
  }
  
};


