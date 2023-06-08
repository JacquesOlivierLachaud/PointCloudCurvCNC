/**
 Copyright (c) 2021
 Jacques-Olivier Lachaud (\c jacques-olivier.lachaud@univ-savoie.fr)
 Laboratory of Mathematics (CNRS, UMR 5127), University of Savoie, France,
 
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
#include <cmath>
#include <algorithm>
#include <stack>
#include <vector>
#include <limits>

/**
 * @file LinearKDTree.h
 * @author Jacques-Olivier Lachaud (\c jacques-olivier.lachaud@univ-savoie.fr)
 * Laboratory of Mathematics (CNRS, UMR 5127), University of Savoie, France
 *
 * @date 2021/04/03
 *
 */

/**
 * @brief A k-d-tree structure stored as an array. The k-d-tree is a
 * strict order that is not computable locally.  This structure is
 * useful for computing nearby points.
 *
 * @tparam TPoint a model for points, i.e. any type with an array subscript operator
 * @tparam dimension the dimension of each point, i.e. its array size.
 */
template < typename TPoint,
           int dimension = TPoint::dimension >
struct LinearKDTree {
  typedef TPoint                   Point;
  typedef std::size_t              Index;
  typedef std::size_t              Size;
  typedef std::vector< Index >     Indices;
  typedef std::vector< Point >     Points;
  typedef double                   Scalar;

  // ------------------------------ public data ------------------------------
public:
  
  Points  _points;  ///< the array of points (stored as given to \ref init).
  Indices _indices; ///< the permutation of \ref _points that gives the order in the tree.

  // ------------------------------ standard services ------------------------------
public:
  /// @name Standard services
  /// @{  
  
  /// Default constructor. Does nothing.
  LinearKDTree() {}

  /// Constructor from a vector of points.
  /// @param points the vector of points that will be structured in the k-d-tree.
  LinearKDTree( const Points& points )
  {
    init( points );
  }

  /// @return the number of points stored in this object.
  Size size() const
  {
    return _points.size();
  }
  
  /// Initializes tge k-d-tree from the given vector of points.
  /// @param points any vector of points.
  void init( const Points& points )
  {
    _points  = points;
    _indices = std::vector<Index>( points.size() );
    for ( Index i = 0; i < points.size(); i++ ) _indices[ i ] = i;
    buildKDTree( 0, points.size(), 0 );
  }

  /// @return the vector of points
  const Points& positions() const
  {
    return _points;
  }

  /// @param v any index of point
  /// @return the corresponding point
  Point position(const Index v) const
  {
    return _points[ v ];
  }

  /// @}
  
  // ------------------------------ localization services ------------------------------
public:
  /// @name Localization services
  /// @{  

  /// This method is useful for guessing a ball radius r, such that when
  /// you look for \a knn neighbors for a point \a p, you have
  /// approximately the same number in the ball of radius r centered on
  /// \a p.
  ///
  /// @param knn the desired number of neighbors.
  ///
  /// @param nb_tests the number of samples that is used to estimate the
  /// radius (if 0, takes the squared root of the number of points).
  ///
  /// @return the estimated radius.
  double findRadius( int knn, int nb_tests = 0 ) const
  {
    RealPoint lo = _points[ _indices.front() ];
    RealPoint up = _points[ _indices.back() ];
    double diameter = ( up - lo ).norm();
    int   number = 0;
    double sum_rho    = 0.0;
    double sum_sqrt_number = 0.0;
    const int nb = _points.size();
    if ( nb_tests <= 0 ) nb_tests = (int) ceil( sqrt( _points.size() ) );
    for ( int i = 0; i < nb_tests; i++ ) {
      RealPoint p = _points[ i % nb ];
      double rho = diameter / sqrt( (double) nb );
      do {
	auto near_points = pointsInBall( p, rho );
	number = near_points.size();
	if ( number >= 3 ) {
	  sum_rho += rho;
	  sum_sqrt_number += sqrt( number );
	}
	rho *= 1.5;
      } while ( number < 30 );
    }
    double best_rho = sqrt( (double) knn ) * sum_rho / sum_sqrt_number;
    // std::cout << "Found rho=" << best_rho << " for " << knn << " neighbors"
    // 	      << ", estimated with " << nb_tests << " samples."
    // 	      << std::endl;
    return best_rho;
  }

  
  /// Given a point \a p, returns the index of its closest neighbor as
  /// well as its squared distance.
  ///
  /// @param[in] p any point in the space
  ///
  /// @return a tuple `(idx,d2)` that is the index of its closest
  /// neighbor and its squared distance to \a p.
  ///
  /// @note This method is put here for completion of localization
  /// services, but remains 6-8x times slower than nanoflann's
  /// findNeighbors.
  std::pair< Index, Scalar >
  nearestNeighbor( const Point p ) const
  {
    typedef std::tuple< Index, Index, int > Node;
    Index  best_im = -1;
    Scalar best_d2 = std::numeric_limits<Scalar>::infinity();
    if ( _indices.empty() ) return std::make_pair( best_im, best_d2 );
    // Looks ugly, but allows searches in vector of points of 2^100 points.
    Node VQ[ 100 ];
    Index  i = 0;
    Index  j = _indices.size();
    int    a = 0;
    Index  m = (i+j)/2;
    Index  im = _indices[ m ];
    best_d2   = distance2( p, _points[ im ] );
    best_im   = im;
    VQ[ 0 ] = { i, j, 0 };
    std::size_t t = 1;
    // Essentially a depth-first binary tree exploration.
    while ( t != 0 ) {
      t -= 1; 
      std::tie( i, j, a ) = VQ[ t ];
      if ( i >= j ) continue;
      // If i and j are close, we switch to a stupid array scan.
      if ( j - i < 8 ) {
        for ( auto k = i; k < j; k++ ) {
          im = _indices[ k ];
          Scalar d2 = distance2( p, _points[ im ] );
          if ( d2 < best_d2 ) {
            best_d2 = d2;
            best_im = im;
          }
        }
      } else {
	// Otherwise it is more a dichotomy.
        m  = (i+j)/2;
        im = _indices[ m ];
        const Scalar dv = _points[ im ][ a ] - p[ a ];
        const auto   na = (a+1) % dimension;
        if ( dv*dv < best_d2 )
          {
            Scalar d2 = distance2( p, _points[ im ] );
            if ( d2 < best_d2 ) {
              best_d2 = d2;
              best_im = im;
            }
            if ( dv > 0.0 ) {
              VQ[ t++ ] = { m+1, j, na };          
              VQ[ t++ ] = { i, m, na };
            } else {
              VQ[ t++ ] = { i, m, na };
              VQ[ t++ ] = { m+1, j, na };          
            }
          }
        else if ( dv > 0.0 )
          VQ[ t++ ] = { i, m, na };
        else
          VQ[ t++ ] = { m+1, j, na };
      }
    }
    return std::make_pair( best_im, best_d2 );
  }

  /// Localization query that returns all the points within a ball.
  ///
  /// @param[in] p any point of the space
  /// @param[in] rho any non-negative value
  /// @param[in] nb_expected a hint for the number of expected points
  /// within the specified ball.
  ///
  /// @return the indices of all the points within the ball of center
  /// \a p and radius \a rho.
  ///
  /// @note as fast as nanoflann's findNeighbors.
  Indices
  pointsInBall( const Point p, const Scalar rho, const Size nb_expected = 50 ) const
  {
    typedef std::tuple< Index, Index, int > Node;
    const Scalar rho2 = rho * rho;
    Indices output;
    output.reserve( nb_expected );
    if ( _indices.empty() ) return output;
    // Looks ugly, but allows searches in vector of points of 2^100 points.
    Node VQ[ 100 ];
    std::size_t t = 1;
    VQ[ 0 ] = { 0, _indices.size(), 0 };
    Index i, j, m, im;
    int   a;
    // Essentially a depth-first binary tree exploration.
    while ( t != 0 ) {
      t -= 1;
      std::tie( i, j, a ) = VQ[ t ];
      if ( i >= j ) continue;
      m  = (i+j)/2;
      im = _indices[ m ];
      if ( distance2( p, _points[ im ] ) <= rho2 )
        output.push_back( im );
      if ( i+1 >= j ) continue;
      const auto na = (a+1) % dimension;
      const auto v  = _points[ im ][ a ];
      if ( ( p[ a ] >= v - rho ) )
        VQ[ t++ ] = { m+1, j, na };
      if ( ( p[ a ] <= v + rho ) )
        VQ[ t++ ] = { i, m, na };
    }
    return output;
  }

  /// Localization query that returns at least the \a k nearest
  /// neighbors. Set \a force_sort to true so that the returned \a k
  /// first are indeed the \a k nearest neighbors.
  ///
  /// @param p any point of the space
  /// @param k the minimum number of returned nearest neighbors
  /// @param rho the starting ball around \a p for searching for
  /// neighbors, a hint that must be positive.
  /// @param force_sort when 'true' the returned vector contains the
  /// points in increasing distance order to \a p, hence the \a k
  /// first are indeed the \a k nearest ones.
  ///
  /// @return the indices of at least \a k points that are the closest to \a p.
  Indices
  kNeighborsAtLeast( const Point& p, int k, Scalar rho, bool force_sort=false) const
  {
    static const Scalar sqrt_2 = sqrt( 2.0 );
    Indices output = pointsInBall( p, rho );
    while ( output.size() < k ) {
      rho   *= sqrt_2;
      output = pointsInBall( p, rho );
    }
    
    // Sorting by distance to retrieve the exact k-nn
    if (force_sort)
      std::sort( output.begin(), output.end(),
		 [&]( const Index& a, const Index& b) {
		   return distance2(p, position(a)) < distance2(p, position(b)); }
		 );
    return output;
  }

  /// @}

  // ---------------------------- basic services ------------------------------
public:
  /// @name Basic services
  /// @{  

  /// @param x any number.
  /// @return its square `x*x`.
  static
  Scalar
  sqr( Scalar x ) const
  {
    return x * x;
  }

  /// @param[in] p,q any two points.
  /// @return their squared Euclidean distance.
  static
  Scalar
  distance2( const Point& p, const Point& q ) const
  {
    Scalar sum = 0.0;
    for ( int i = 0; i < dimension; i++ )
      sum += sqr( p[ i ] - q[ i ] );
    return sum;
  }

  /// @}

  // ---------------------------- internal methods ------------------------------
protected:

  // Internal method that builds the k-d-tree recursively.
  void buildKDTree( Index i, Index j, int a )
  {
    if ( i+1 >= j ) return;
    Index m = (i+j)/2;
    std::nth_element( _indices.begin() + i, _indices.begin() + m, _indices.begin() + j,
                      [&] ( Index ip, Index jp )
                      { return _points[ ip ][ a ] < _points[ jp ][ a ]; } );
    buildKDTree( i,   m, (a+1) % dimension );
    buildKDTree( m+1, j, (a+1) % dimension );
  }
  
};
