#include "polyscope/polyscope.h"
#include "polyscope/point_cloud.h"
#include "Common.h"
#include "CorrectedNormalCurrentFormulaEigen.h"
#include "LinearKDTree.h"
#include "PointCloudCurvatureComputer.h"
#include "RegularPolyhedron.h"
#include "Helpers.h"

typedef Eigen::VectorXd  DenseVector;
typedef Eigen::MatrixXd  DenseMatrix;
typedef Eigen::Vector3d  RealPoint;
typedef Eigen::Vector3d  RealVector;
typedef Eigen::Matrix3d  RealTensor;
typedef PointCloudCurvatureComputer::Method GenerationMethod;

// Polyscope visualization handle, to quickly add data to the surface
polyscope::PointCloud* ptCloud = nullptr;
polyscope::PointCloudScalarQuantity* ptCloudH = nullptr;
polyscope::PointCloudScalarQuantity* ptCloudG = nullptr;
polyscope::PointCloudScalarQuantity* ptCloudK1 = nullptr;
polyscope::PointCloudScalarQuantity* ptCloudK2 = nullptr;
// Global variables for user interaction
float Time         = 1.0;
double TotalTime   = 0.0;
int number_of_points = 10000;
int   nn_normals   = 20;
int   k_nn  = 20;
int   nn_density   = 20;
int   maxtriangles = 200;
float sdistance    = 0.1;
float epsilon      = 0.0;
float xi           = 0.0;
int   scale        = 4;
float AvgNormalsWeight  = 0.5f;
bool FastRenderMode     = true;
float radius_R     = 1.0;
float radius_r     = 0.3;

// Main global variables
PointCloudCurvatureComputer  cnc_computer;
std::string                  fname;
std::vector< RealPoint >     points;
std::vector< RealVector >    point_normals;
std::vector< RealPoint >     obj_points;
std::vector< RealVector >    obj_normals;
LinearKDTree< RealPoint, 3 > ptTree;
GenerationMethod             method;

// bounding box
RealPoint lowest;
RealPoint uppest;

std::vector< double > true_H;
std::vector< double > true_G;
std::vector< double > true_K1;
std::vector< double > true_K2;

double H_error_l2   = 0.0;
double H_error_loo  = 0.0;
double G_error_l2   = 0.0;
double G_error_loo  = 0.0;
double K1_error_l2  = 0.0;
double K1_error_loo = 0.0;
double K2_error_l2  = 0.0;
double K2_error_loo = 0.0;

// -------------------- ERRORS --------------------
template < typename A1, typename A2 >
double error_l2( int n, const A1& a1, const A2& a2 )
{
  double sum = 0.0;
  for ( int i = 0; i < n; i++ )
    sum += ( a1[ i ] - a2[ i ] ) * ( a1[ i ] - a2[ i ] );
  return sqrt( sum / n );
}
template < typename A1, typename A2 >
double error_loo( int n, const A1& a1, const A2& a2 )
{
  double m = 0.0;
  for ( int i = 0; i < n; i++ )
    m = std::max( m, fabs( a1[ i ] - a2[ i ] ) );
  return m;
}

std::pair< RealPoint, RealPoint >
computeBoundingBox( const std::vector< RealPoint >& pts )
{
  RealPoint lo = pts[ 0 ];
  RealPoint hi = pts[ 0 ];
  for ( auto&& p : pts ) {
    for ( int i = 0; i < 3; i++ ) {
      lo[ i ] = std::min( lo[ i ], p[ i ] );
      hi[ i ] = std::max( hi[ i ], p[ i ] );
    }
  }
  return std::make_pair( lo, hi );
}


// -------------------- CURVATURES --------------------                     
void displayCurvatures
( const DenseVector& mA,   const DenseVector& mH,   const DenseVector& mG,
  const DenseVector& mT11, const DenseVector& mT12, const DenseVector& mT13,
  const DenseVector& mT22, const DenseVector& mT23, const DenseVector& mT33 )
{
  const auto nb = mA.rows();
  DenseVector H( nb );
  DenseVector G( nb );
  DenseVector K1( nb );
  DenseVector K2( nb );
  std::vector< RealVector > V1( nb );
  std::vector< RealVector > V2( nb );
  double     k1, k2;
  RealVector v1, v2;
  for ( int i = 0; i < nb; i++ )
    {
      const RealPoint p = points[ i ];
      H[ i ] = ( mA[ i ] == 0.0 ) ? 0.0 : ( mH[ i ] / mA[ i ] );
      G[ i ] = ( mA[ i ] == 0.0 ) ? 0.0 : ( mG[ i ] / mA[ i ] );
      RealTensor T;
      if ( mA[ i ] != 0.0 )
	{
	  T <<
	    mT11[ i ], mT12[ i ], mT13[ i ],
	    mT12[ i ], mT22[ i ], mT23[ i ],
	    mT13[ i ], mT23[ i ], mT33[ i ];
	  T /= mA[ i ];
	}
      std::tie( K1[ i ], K2[ i ], V1[ i ], V2[ i ] ) = 
        CNCEigen::curvaturesFromTensor( T, 1.0, point_normals[ i ] );
    }
  if ( ptCloudH != nullptr )
    {
      const auto range = ptCloudH->getMapRange();
      ptCloudH = ptCloud->addScalarQuantity( "Curvature H (P)", H );
      ptCloudH->setMapRange( range );
    }
  else ptCloudH = ptCloud->addScalarQuantity( "Curvature H (P)", H );
  if ( ptCloudG != nullptr )
    {
      const auto range = ptCloudG->getMapRange();
      ptCloudG = ptCloud->addScalarQuantity( "Curvature G (P)", G );
      ptCloudG->setMapRange( range );
    }
  else ptCloudG = ptCloud->addScalarQuantity( "Curvature G (P)", G );
  if ( ptCloudK1 != nullptr )
    {
      const auto range = ptCloudK1->getMapRange();
      ptCloudK1 = ptCloud->addScalarQuantity( "Curvature K1 (P)", K1 );
      ptCloudK1->setMapRange( range );
    }
  else ptCloudK1 = ptCloud->addScalarQuantity( "Curvature K1 (P)", K1 );
  if ( ptCloudK2 != nullptr )
    {
      const auto range = ptCloudK2->getMapRange();
      ptCloudK2 = ptCloud->addScalarQuantity( "Curvature K2 (P)", K2 );
      ptCloudK2->setMapRange( range );
    }
  else ptCloudK2 = ptCloud->addScalarQuantity( "Curvature K2 (P)", K2 );
  ptCloud->addVectorQuantity( "Curvature V1 (P)", V1 );
  ptCloud->addVectorQuantity( "Curvature V2 (P)", V2 );

  if ( true_H.size() == H.rows() )
    {
      std::cout << "Computing H errors" << std::endl;
      H_error_l2  = error_l2 ( H.rows(), true_H, H );
      H_error_loo = error_loo( H.rows(), true_H, H );
    }
  if ( true_G.size() == G.rows() )
    {
      G_error_l2  = error_l2 ( G.rows(), true_G, G );
      G_error_loo = error_loo( G.rows(), true_G, G );
    }
  if ( true_K1.size() == K1.rows() )
    {
      K1_error_l2  = error_l2 ( K1.rows(), true_K1, K1 );
      K1_error_loo = error_loo( K1.rows(), true_K1, K1 );
    }
  if ( true_K2.size() == K2.rows() )
    {
      K2_error_l2  = error_l2 ( K2.rows(), true_K2, K2 );
      K2_error_loo = error_loo( K2.rows(), true_K2, K2 );
    }

  
}


void doCurvatures()
{
  if ( point_normals.size() != ptTree.size() ) return;
  std::cout << "Computing curvatures of " << ptTree.size() << " points." << std::endl;
  cnc_computer.init( ptTree, k_nn );
  if ( method == GenerationMethod::IndependentGeneration )
    cnc_computer.chooseIndependentGenerationMethod( maxtriangles );
  else if ( method == GenerationMethod::HexagramGeneration )
    cnc_computer.chooseHexagramGenerationMethod( AvgNormalsWeight );
  else if ( method == GenerationMethod::AvgHexagramGeneration )
    cnc_computer.chooseAvgHexagramGenerationMethod( AvgNormalsWeight );
  else
    cnc_computer.chooseUniformGenerationMethod( maxtriangles );
  
  DenseVector A, H, G, T11, T12, T13, T22, T23, T33;
  std::vector<double> timings_cnc;
  std::tie( A, H, G, T11, T12, T13, T22, T23, T33, timings_cnc) =
    cnc_computer.computeCurvatureMeasures( point_normals );
  TotalTime = 0.0;
  for ( auto t : timings_cnc ) TotalTime += t;
  TotalTime /= 1e6; // time in ms
  displayCurvatures( A, H, G, T11, T12, T13, T22, T23, T33 );
}

// -------------------- NORMALS --------------------
// void doTrueNormals()
// {
//   // compute normals
//   if ( true_normals.size() != ptTree._points.size() ) return;
//   double rho = findRadius( ptTree, nn_normals );
//   // barycenter_normals = computeBarycenterNormals( nn_normals, rho );
//   // ptCloud->addVectorQuantity( "Barycenter normals", barycenter_normals );
//   point_normals = true_normals;
//   ptCloud->addVectorQuantity( "Point normals", point_normals );
// }


// -------------------- SHAPE GENERATION --------------------
void clearTrueValues()
{
  true_H.clear();
  true_G.clear();
  true_K1.clear();
  true_K2.clear();
  H_error_l2  = 0.0;
  H_error_loo = 0.0;
  G_error_l2  = 0.0;
  G_error_loo = 0.0;
  K1_error_l2  = 0.0;
  K1_error_loo = 0.0;
  K2_error_l2  = 0.0;
  K2_error_loo = 0.0;
}

std::pair< std::vector< RealPoint >, std::vector< RealVector > >
generateCube( int N, float R = 1.0 )
{
  std::vector< RealPoint >  points;
  std::vector< RealVector > true_normals;
  std::cout << "Build " << N << " points ..." << std::endl;
  RegularPolyhedron< RealPoint > cube( 6, R );
  for (size_t i = 0; i < N; i++) {
    RealPoint p;
    RealVector n;
    std::tie( p, n ) = cube.randomPointNormal();
    points.push_back( p );
    true_normals.push_back( n );
  }
  clearTrueValues();
  return std::make_pair( points, true_normals );
}

std::pair< std::vector< RealPoint >, std::vector< RealVector > >
generateDodecahedron( int N, float R = 1.0 )
{
  std::vector< RealPoint >  points;
  std::vector< RealVector > true_normals;
  std::cout << "Build " << N << " points ..." << std::endl;
  RegularPolyhedron< RealPoint > cube( 12, R );
  for (size_t i = 0; i < N; i++) {
    RealPoint p;
    RealVector n;
    std::tie( p, n ) = cube.randomPointNormal();
    points.push_back( p );
    true_normals.push_back( n );
  }
  clearTrueValues();
  return std::make_pair( points, true_normals );
}

std::pair< std::vector< RealPoint >, std::vector< RealVector > >
generateSphere( int N, float R = 1.0 )
{
  std::vector< RealPoint >  points;
  std::vector< RealVector > true_normals;
  std::cout << "Build " << N << " points ..." << std::endl;
  true_H.resize( N );
  true_G.resize( N );
  true_K1.resize( N );
  true_K2.resize( N );
  for (size_t i = 0; i < N; i++) {
    glm::vec3 u = { polyscope::randomUnit() - .5, 
      polyscope::randomUnit() - .5, 
      polyscope::randomUnit() - .5 };
    glm::vec3 n = u / glm::length( u );
    u = ( R / glm::length( u ) ) * u;
    points.push_back( RealPoint( u.x, u.y, u.z ) );
    true_normals.push_back( RealVector( n.x, n.y, n.z ) );
    true_H[ i ] = 1.0 / R;
    true_G[ i ] = 1.0 / (R*R);
    true_K1[ i ] = 1.0 / R;
    true_K2[ i ] = 1.0 / R;
  }
  return std::make_pair( points, true_normals );
}

std::pair< std::vector< RealPoint >, std::vector< RealVector > >
generateTorus( int N, float R = 1.0, float r = 0.2 )
{
  std::vector< RealPoint > points;
  std::vector< RealVector > true_normals;
  std::cout << "Build " << N << " points ..." << std::endl;
  true_H.resize( N );
  true_G.resize( N );
  true_K1.resize( N );
  true_K2.resize( N );
  for (size_t i = 0; i < N; i++) {
    glm::vec3 u = { polyscope::randomUnit() - .5, polyscope::randomUnit() - .5, 0 };
    glm::vec3 v = { 0.0, 0.0, 1.0 };
    u = u / glm::length( u ); // direction
    float theta = polyscope::randomUnit() * 2.0 * M_PI;
    glm::vec3 b = R * u;
    glm::vec3 x = b + r * ( cos(theta) * u + sin(theta) * v );
    glm::vec3 n = cos(theta) * u + sin(theta) * v;
    points.push_back( RealPoint( x.x, x.y, x.z ) );
    true_normals.push_back( RealVector( n.x, n.y, n.z ) );
    true_K1[ i ] = 1.0 / r;
    true_K2[ i ] = cos( theta ) / (R+r*cos(theta));
    true_H[ i ] = 0.5 * ( true_K1[ i ] + true_K2[ i ] );
    true_G[ i ] = true_K1[ i ] * true_K2[ i ];
  }
  return std::make_pair( points, true_normals );
}

std::vector< RealPoint >
hausdorffPerturbation( const std::vector< RealPoint >& pts, double eps )
{
  std::vector< RealPoint > output = pts;
  glm::vec3 u;
  for ( auto& p : output ) {
    do {
      u = { polyscope::randomUnit() - .5,
        polyscope::randomUnit() - .5, polyscope::randomUnit() - .5 };
    } while ( glm::length( u ) >= 0.5 );
    u *= (float) (2.0*eps);
    p = RealPoint( p[ 0 ]+u.x, p[ 1 ]+u.y, p[ 2 ]+u.z );
  }
  return output;
}

std::vector< RealPoint >
uniformNormalPerturbation( const std::vector< RealVector >& n, double eps )
{
  std::vector< RealVector > output = n;
  glm::vec3 u;
  for ( auto& p : output ) {
    do {
      u = { polyscope::randomUnit() - .5,
        polyscope::randomUnit() - .5, polyscope::randomUnit() - .5 };
    } while ( glm::length( u ) >= 0.5 );
    u *= (float) (2.0*eps);
    p = RealVector( p[ 0 ]+u.x, p[ 1 ]+u.y, p[ 2 ]+u.z );
    p /= p.norm();
  }
  return output;
}
    
void doGenerateCube()
{
  ptCloudH = nullptr;
  ptCloudG = nullptr;
  ptCloudK1 = nullptr;
  ptCloudK2 = nullptr;
  // generate points
  std::tie( points, point_normals )
    = generateCube( number_of_points, radius_R );
  points = hausdorffPerturbation( points, epsilon );
  point_normals = uniformNormalPerturbation( point_normals, xi );
  if ( ptCloud != nullptr )
    polyscope::removePointCloud( "point cloud" );
  std::cout << "Build k-d-tree ..." << std::endl;
  ptTree.init( points );
  // visualize!
  std::cout << "register point cloud ..." << std::endl;
  ptCloud = polyscope::registerPointCloud("point cloud", points);
  ptCloud->setPointRadius( radius_R / 100.0, false );
  ptCloud->setPointRenderMode( FastRenderMode
                               ? polyscope::PointRenderMode::Quad
                               : polyscope::PointRenderMode::Sphere );
  std::cout << "... done" << std::endl;
  std::tie( lowest, uppest ) = computeBoundingBox( points );
}

void doGenerateDodecahedron()
{
  ptCloudH = nullptr;
  ptCloudG = nullptr;
  ptCloudK1 = nullptr;
  ptCloudK2 = nullptr;
  // generate points
  std::tie( points, point_normals )
    = generateDodecahedron( number_of_points, radius_R );
  points = hausdorffPerturbation( points, epsilon );
  point_normals = uniformNormalPerturbation( point_normals, xi );
  if ( ptCloud != nullptr )
    polyscope::removePointCloud( "point cloud" );
  std::cout << "Build k-d-tree ..." << std::endl;
  ptTree.init( points );
  // visualize!
  std::cout << "register point cloud ..." << std::endl;
  ptCloud = polyscope::registerPointCloud("point cloud", points);
  ptCloud->setPointRadius( radius_R / 100.0, false );
  ptCloud->setPointRenderMode( FastRenderMode
                               ? polyscope::PointRenderMode::Quad
                               : polyscope::PointRenderMode::Sphere );
  std::cout << "... done" << std::endl;
  std::tie( lowest, uppest ) = computeBoundingBox( points );
}

void doGenerateSphere()
{
  ptCloudH = nullptr;
  ptCloudG = nullptr;
  ptCloudK1 = nullptr;
  ptCloudK2 = nullptr;
  // generate points
  std::tie( points, point_normals )
    = generateSphere( number_of_points, radius_R );
  points = hausdorffPerturbation( points, epsilon );
  point_normals = uniformNormalPerturbation( point_normals, xi );
  if ( ptCloud != nullptr )
    polyscope::removePointCloud( "point cloud" );
  std::cout << "Build k-d-tree ..." << std::endl;
  ptTree.init( points );
  // visualize!
  std::cout << "register point cloud ..." << std::endl;
  ptCloud = polyscope::registerPointCloud("point cloud", points);
  ptCloud->setPointRadius( radius_R / 100.0, false );
  ptCloud->setPointRenderMode( FastRenderMode
                               ? polyscope::PointRenderMode::Quad
                               : polyscope::PointRenderMode::Sphere );
  std::cout << "... done" << std::endl;
  std::tie( lowest, uppest ) = computeBoundingBox( points );
}

void doGenerateTorus()
{
  ptCloudH = nullptr;
  ptCloudG = nullptr;
  ptCloudK1 = nullptr;
  ptCloudK2 = nullptr;

  // generate points
  std::tie( points, point_normals )
    = generateTorus( number_of_points, radius_R, radius_r );
  points = hausdorffPerturbation( points, epsilon );
  point_normals = uniformNormalPerturbation( point_normals, xi );
  if ( ptCloud != nullptr )
    polyscope::removePointCloud( "point cloud" );
  std::cout << "Build k-d-tree ..." << std::endl;
  ptTree.init( points );
  // visualize!
  std::cout << "register point cloud ..." << std::endl;
  ptCloud = polyscope::registerPointCloud("point cloud", points);
  ptCloud->setPointRadius( radius_R / 100.0, false );
  ptCloud->setPointRenderMode( FastRenderMode
                               ? polyscope::PointRenderMode::Quad
                               : polyscope::PointRenderMode::Sphere );
  std::cout << "... done" << std::endl;
  std::tie( lowest, uppest ) = computeBoundingBox( points );
}

void doGenerateFromFile()
{
  ptCloudH  = nullptr;
  ptCloudG  = nullptr;
  ptCloudK1 = nullptr;
  ptCloudK2 = nullptr;
  if ( fname == "" ) return;
  std::ifstream input( fname );
  std::string extname = fname.substr( fname.find_last_of(".") + 1 );
  std::tie( obj_points, obj_normals )
    = extname == "obj"
    ? readOBJAsPointCloud<RealPoint,RealVector>( input )
    : readPointCloud<RealPoint,RealVector>     ( input );
  input.close();
  points = hausdorffPerturbation( obj_points, epsilon );
  point_normals = uniformNormalPerturbation( obj_normals, xi );  
  if ( ptCloud != nullptr )
    polyscope::removePointCloud( "point cloud" );
  std::cout << "Build k-d-tree ..." << std::endl;
  ptTree.init( points );
  // visualize!
  std::cout << "register point cloud ..." << std::endl;
  ptCloud = polyscope::registerPointCloud("point cloud", points);
  ptCloud->setPointRadius( radius_R / 100.0, false );
  ptCloud->setPointRenderMode( FastRenderMode
                               ? polyscope::PointRenderMode::Quad
                               : polyscope::PointRenderMode::Sphere );
  std::cout << "... done" << std::endl;
  std::tie( lowest, uppest ) = computeBoundingBox( points );
}



// -------------------- IHM --------------------                     
void myCallback()
{
  ImGui::Text("Generated input data");
  ImGui::SameLine();
  if ( ImGui::Checkbox("Fast Render", &FastRenderMode ) )
    ptCloud->setPointRenderMode( FastRenderMode
				 ? polyscope::PointRenderMode::Quad
				 : polyscope::PointRenderMode::Sphere );

  ImGui::SliderInt("N (nb of samples)", &number_of_points, 1000, 1000000 );
  ImGui::SliderFloat("R (big radius)", &radius_R, 0.0, 1.0);
  ImGui::SliderFloat("r (small radius)", &radius_r, 0.0, 1.0);
  ImGui::SliderFloat("x perturbation", &epsilon, 0.0, 0.2);
  ImGui::SliderFloat("u perturbation", &xi, 0.0, 0.2);
  if (ImGui::Button("Sphere")) doGenerateSphere();
  ImGui::SameLine();
  if (ImGui::Button("Torus"))  doGenerateTorus();
  ImGui::SameLine();
  if (ImGui::Button("Cube"))   doGenerateCube();
  ImGui::SameLine();
  if (ImGui::Button("Dodecahedron")) doGenerateDodecahedron();
  ImGui::SameLine();
  if (ImGui::Button("InputFile"))    doGenerateFromFile();
  ImGui::Text("CNC Curvatures estimation");
  ImGui::SliderInt("K (#nearest neighbors)", &k_nn, 3, 100);
  ImGui::SliderInt("T (max triangles)", &maxtriangles, 1, 1000);
  if ( ImGui::RadioButton( "Uniform",
                           method == GenerationMethod::UniformGeneration) )
    method = GenerationMethod::UniformGeneration;
  ImGui::SameLine();
  if ( ImGui::RadioButton( "Independent",
                           method == GenerationMethod::IndependentGeneration ) )
    method = GenerationMethod::IndependentGeneration;
  ImGui::SameLine();
  if ( ImGui::RadioButton( "Hexagram",
                           method == GenerationMethod::HexagramGeneration ) )
    method = GenerationMethod::HexagramGeneration;
  ImGui::SameLine();
  if ( ImGui::RadioButton( "Avg-Hexagram",
                           method == GenerationMethod::AvgHexagramGeneration ) )
    method = GenerationMethod::AvgHexagramGeneration;
  ImGui::SliderFloat("Avg normals weight", &AvgNormalsWeight, 0.0f, 1.0f);
  if (ImGui::Button("Curvatures")) doCurvatures(); // either curvatures or octree curvatures
  ImGui::SameLine();
  ImGui::Text( "Total time=%lf ms", TotalTime );

  ImGui::Text( "H l2 error = %f, H loo error = %f",
               (double) H_error_l2, (double) H_error_loo );
  ImGui::Text( "G l2 error = %f, G loo error = %f",
               (double) G_error_l2, (double) G_error_loo );
  ImGui::Text( "K1 l2 error = %f, K1 loo error = %f",
               (double) K1_error_l2, (double) K1_error_loo );
  ImGui::Text( "K2 l2 error = %f, K2 loo error = %f",
               (double) K2_error_l2, (double) K2_error_loo );
  
}

int main(int argc, char **argv)
{
  polyscope::init();

  // Set the callback function
  polyscope::state::userCallback = myCallback;

  // Process input file if any
  fname = argc > 1 ? argv[ 1 ] : "";
  if ( fname != "" )
    doGenerateFromFile();
  else
    doGenerateTorus();

  method = GenerationMethod::AvgHexagramGeneration;
  
  // Run polyscope viewer
  polyscope::show();
  
  return 0;
}

