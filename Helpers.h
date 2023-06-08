#pragma once

#include "Common.h"


//-----------------------------------------------------------------------------
template <typename RealPoint, typename RealVector>
std::pair< std::vector<RealPoint>, std::vector<RealVector> >
readOBJAsPointCloud( std::istream & input )
{
  std::vector<RealPoint>  vertices;
  std::vector<RealVector> normals;
  std::string linestr;
  std::string keyword;
  std::string indices;
  RealPoint  p;
  RealVector n;
  std::getline( input, linestr );
  int l = 0;
  for ( ; input.good() && ! input.eof(); std::getline( input, linestr ), l++ )
  {
    if ( linestr.empty() ) continue; // skip empty line
    if ( linestr[0] == '#' ) continue; // skip comment line
    std::istringstream lineinput( linestr );
    std::operator>>( lineinput, keyword ); // lineinput >> keyword;
    if ( keyword == "v" ) {
      lineinput >> p[ 0 ] >> p[ 1 ] >> p[ 2 ];
      vertices.push_back( p );
    } else if ( keyword == "vn" ) {
      lineinput >> n[ 0 ] >> n[ 1 ] >> n[ 2 ];
      normals.push_back( n );
    }
    // Weird: necessary to clear them.
    keyword = ""; linestr = "";
  }
  // Creating SurfaceMesh
  std::cout << "[readOBJAsPointCloud] Read"
  << " #lines=" << l
  << " #V=" << vertices.size()
  << " #VN=" << normals.size()
  << std::endl;
  if ( input.bad() )
    std::cerr << "[readOBJAsPointCloud] Some I/O error occured."
    << " Proceeding but the point cloud  may be damaged." << std::endl;
  return std::make_pair( vertices, normals );
}
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Read a file 
template <typename RealPoint, typename RealVector>
std::pair< std::vector<RealPoint>, std::vector<RealVector> >
readPointCloud( std::istream & input )
{
  std::vector<RealPoint>  vertices;
  std::vector<RealVector> normals;
  std::string linestr;
  std::string indices;
  RealPoint  p;
  RealVector n;
  std::getline( input, linestr );
  int l = 0;
  for ( ; input.good() && ! input.eof(); std::getline( input, linestr ), l++ )
  {
    if ( linestr.empty() ) continue; // skip empty line
    if ( linestr[0] == '#' ) continue; // skip comment line
    std::istringstream lineinput( linestr );
    lineinput >> p[ 0 ] >> p[ 1 ] >> p[ 2 ];
    vertices.push_back( p );
    lineinput >> n[ 0 ] >> n[ 1 ] >> n[ 2 ];
    normals.push_back( n );
    linestr = "";
  }
  // Creating SurfaceMesh
  std::cout << "[readPointCloud] Read"
  << " #lines=" << l
  << " #V=" << vertices.size()
  << " #VN=" << normals.size()
  << std::endl;
  if ( input.bad() )
    std::cerr << "[readPointCloud] Some I/O error occured."
	      << " Proceeding but the point cloud  may be damaged." << std::endl;
  return std::make_pair( vertices, normals );
}
//-----------------------------------------------------------------------------


//-----------------------------------------------------------------------------
template <typename RealPoint>
std::vector<RealPoint> addUniformNoise( const std::vector<RealPoint> &input,
                                       double eps)
{
  //hardcoded noise
  std::random_device rd;  //Will be used to obtain a seed for the random number engine
  std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
  std::uniform_real_distribution<> dis(0.0, 1.0);
  std::vector<RealPoint> points(input.size());

  auto cpt=0u;
  for(auto &p: input)
  {
    RealPoint dec = {eps*dis(gen),eps*dis(gen),eps*dis(gen)};
    points[cpt] = p + dec;
    cpt++;
  }
  return points;
}
//-----------------------------------------------------------------------------
