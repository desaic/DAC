#ifndef ELEMENTMESH_HPP
#define ELEMENTMESH_HPP

#include <vector>
#include <Eigen/Dense>
#include <Eigen/Sparse>

class Element;
class Material;
class Quadrature;
class TrigMesh;

class ElementMesh{

public:

  ///@brief elements will be deleted by destructor. Do not put same pointer in different meshes.
  ///Make copies if needed.
  std::vector<Element*>e;

  ///@brief deformed positions in t-1. used by dynamics(unit m)
//  std::vector<Vector3f>x0;

  ///@brief deformed positions (unit m)
  std::vector<Eigen::Vector3d>x;
  
  ///@brief vertex positions at rest pose.
  std::vector<Eigen::Vector3d>X;

  ///@brief velocity
  std::vector<Eigen::Vector3d>v;

  ///unprojected velocity
  std::vector<Eigen::Vector3d>vu;

  ///@brief stiffness blocks.
  std::vector< std::map< int, Eigen::MatrixXd > > sblocks;
  ///@brief stiffness pattern as long as topology of the fem mesh does not changes.
  Eigen::SparseMatrix<double> Kpattern;
  ///@brief dimension of the mesh. either 2 or 3
  ///in 2D, the third coordinate is ignored/wasted.
  int dim;

  ///@brief gravitational constant
  Eigen::Vector3d G;

  ///@brief damping constants for mass and stiffness
  float dampAlpha, dampBeta;

  ///@brief special collision vert
  int specialVert;

  std::vector<int> topVerts;

  ///@brief lumped mass
  std::vector<float>mass;
  Eigen::SparseMatrix<double> M;
  float totalMass;
  float totalVol ;

  ///@brief pointers are not freed by destructor.
  std::vector<Material*>m;
  
  ///@brief Color for different materials used in rendering.
  std::vector<Eigen::Vector3f> materialColor;

  ///@brief material for each element
  std::vector<int>me;
  
  ///@brief per element depth
  std::vector<double> depth;

  ///@brief external forces applied to each dof.
  std::vector<Eigen::Vector3d>fe;
  
  ///@brief Contact/Frictional forces at each dof.
  std::vector<Eigen::Vector3d>fc;

  ///@brief scale for drawing forces
  float forceDrawingScale;
  ///@brief scale depth to use as color for rendering.
  float depthColorScale;

  ///@brief fixed degree of freedom. Expected size #verts x dim.
  std::vector<int> fixedDof;

  ///@brief lower bounds for vertex positions
  std::vector<Eigen::Vector3d> lb;

  ///@brief optional mesh for drawing
  std::vector<std::vector<Eigen::Vector3d> > * u;

  ///@brief default constructor
  ElementMesh();
  
  ///@brief load a plain text hex fem mesh
  void load(std::istream & in, float scale=1.0f);

  ///@brief saves topology and the array x in deformed state.
  void saveMesh(std::ostream & out);

  ///@brief utility. call after initializing or changing X and e 
  ///X is copied to x;
  void initArrays();

  void initElements(const Quadrature *q);

  ///@brief Place such that lower left corner is origin.
  void place(Eigen::Vector3d origin);

  ///@brief add a material to the list of materials in this mesh
  void addMaterial(Material*_m);

  ///@brief for debug, check the size of members.
  int check();

  ///@brief place elements on a 2D grid assuming a regular quadrilateral element grid.
  ///@brief grid[ii][jj] stores element index at that grid index. -1 if empty.
  void grid2D(std::vector<std::vector< int> > & grid) const;
  void grid3D(std::vector<std::vector< std::vector<int> > > & grid)const;

  double getAllEnergy(bool addgravity=true);
  double getEnergy(bool addgravity=true);
  
  ///@brief internal elastic energy.
  double getEnergyElastic();

  double getKineticEnergy(bool addgravity=false);
  ///@brief get elastic energy of on element
  double getEnergyEle(int eIdx);

  ///@brief computes internal forces only
  std::vector<Eigen::Vector3d> getForceEle(int eIdx);
  std::vector<Eigen::Vector3d> getForce(bool add_gravity=true);

  Eigen::MatrixXd getStiffness();
  Eigen::MatrixXd getStiffness(int eIdx);

  ///@brief computes elastic force differential i.e. K*dx
  Eigen::VectorXd getdForceEle(int eIdx, const Eigen::VectorXd & dx);
  Eigen::VectorXd getdForce(const Eigen::VectorXd & dx);

  ///@param I row offsets. I.size() = matrix size + 1. I[size()-1]=number of non-zeros.
  ///@param J column indices.
  void stiffnessPattern(std::vector<int> & I, std::vector<int> & J);

  Eigen::SparseMatrix<double> getStiffnessSparse();
  Eigen::SparseMatrix<double> getStiffnessDamping(Eigen::SparseMatrix<double> & D,
    std::vector<float> dcoef);

  ///@brief precompute lumped maasses.
  void computeMass();

  Eigen::Vector3d centerOfMass() const;

  double getTopAng();

  Eigen::Vector3d linearMomentum() const;

  Eigen::Vector3d linearVelocity() const;

  Eigen::Vector3d
  angularVelocity(Eigen::Vector3d center) const;

  Eigen::Matrix3d
  inertiaTensor(Eigen::Vector3d center)const;

  Eigen::Vector3d
  angularMomentum(Eigen::Vector3d center) const;

  virtual ~ElementMesh();

};

void fixRotation(Eigen::SparseMatrix<double> & K, ElementMesh * mesh);
void fixTranslation(Eigen::SparseMatrix<double> & K, ElementMesh * mesh);

#endif
