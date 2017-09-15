#include "SampleForces.hpp"
#include "ArrayUtil.hpp"
#include "Element.hpp"
#include "ElementHex.hpp"
#include "ElementMesh.hpp"
#include "ElementRegGrid.hpp"
#include "ElementMeshUtil.hpp"
#include "femError.hpp"
#include "FileUtil.hpp"
#include "GenerateForces.hpp"
#include "LinSolver.hpp"
#include "LinPardiso.hpp"
#include "Quadrature.hpp"
#include "Material.hpp"
#include "MaterialQuad.hpp"
#include "StepperStatic.hpp"
#include "StrainEneNeo.hpp"
#include "StrainEneANeo.hpp"
#include "World.hpp"
#include <map>

void loadSamplingWorld(const ConfigFile & conf, World * world)
{
  int dim = 3;
  conf.getInt("dim", dim);
  int gridres = 4;
  conf.getInt("gridres", gridres);
  //ElementMesh * em = new ElementRegGrid(gridres, gridres, gridres);
  //scaleMesh(em, 0.1);
  ElementMesh * em = new ElementMesh();
  em->dim = dim; 
  std::string meshFile = conf.getString("meshfile");
  FileUtilIn in(conf.dir + "/" + meshFile);
  em->load(in.in);

  //make materials
  std::string materialFile = "../../models/material3d.txt";
  if (conf.hasOpt("materialfile")){
    materialFile = conf.dir + conf.getString("materialfile");
  }

  FileUtilIn matIn(materialFile);
  if (!matIn.good()){
    return;
  }
  std::vector<double> densities;
  std::vector<StrainEne*> ene = loadMaterials(matIn.in, densities);
  matIn.close();
  const Quadrature* quadrature = &(Quadrature::Gauss2);
  if (dim == 1){
    quadrature = &(Quadrature::Gauss2_2D);
  }
  std::vector<Material *> material(ene.size());
  for (unsigned int ii = 0; ii<material.size(); ii++){
    MaterialQuad * matq = new MaterialQuad();
    if (dim == 2){
      matq->q = quadrature;
      matq->e.resize(matq->q->x.size());
    }
    for (unsigned int jj = 0; jj < matq->e.size(); jj++){
      matq->e[jj] = ene[ii];
    }
    material[ii] = matq;
    material[ii]->density = densities[ii];
    em->addMaterial(matq);
  }
  em->check();

  world->quadrature = quadrature;
  em->initElements(quadrature);
  world->em_.push_back(em);

  //initialize linear solver stiffness matrix pattern
  //stiffness matrix is fixed for all force samples since the mesh is fixed.
  std::vector<int> I, J;
  for (size_t i = 0; i < world->em_.size(); i++){
    ElementMesh * em = world->em_[i];
    if (world->solvers.size() <= i){
      LinPardiso<double>  * solver = new LinPardiso<double>();
      solver->init();
      em->stiffnessPattern(I, J);
      Eigen::SparseMatrix<double> K = em->Kpattern;
      fixTranslation(K, em);
      fixRotation(K, em);
      solver->init(K);
      world->solvers.push_back(solver);
    }
  }
}

//copy the same force of an element to all elements.
//overlapping vertices get multiple copies of the same force from adjacent elements.
void copyForce(ElementMesh * em, const std::vector<float> & force)
{
  for (int i = 0; i < em->fe.size(); i++){
    em->fe[i] = Eigen::Vector3d::Zero();
  }
  int nV = em->e[0]->nV();
  int dim = 3;
  std::vector<Eigen::Vector3d> forceList(nV);
  for (int i = 0; i < nV; i++){
    for (int j = 0; j < dim; j++){
      forceList[i][j] = force[i*dim + j];
    }
  }
  for (size_t i = 0; i < em->e.size(); i++){
    Element * ele = em->e[i];
    for (int j = 0; j < ele->nV(); j++){
      int vidx = ele->at(j);
      em->fe[vidx] += forceList[j];
    }
  }
}

//copy the same force of an element to all elements.
//overlapping vertices get multiple copies of the same force from adjacent elements.
void copyForce(ElementMesh * em, const std::vector<Eigen::Vector3d> & force)
{
  for (int i = 0; i < em->fe.size(); i++){
    em->fe[i] = Eigen::Vector3d::Zero();
  }
  int nV = em->e[0]->nV();
  for (size_t i = 0; i < em->e.size(); i++){
    Element * ele = em->e[i];
    for (int j = 0; j < ele->nV(); j++){
      int vidx = ele->at(j);
      em->fe[vidx] += force[j];
    }
  }
}

double potentialEnergy(ElementMesh * em)
{
  double E = em->getEnergy(false);
  return E;
}

std::vector<double> staticGradient(ElementMesh * em)
{
  std::vector<Eigen::Vector3d> F = em->getForce(false);
  int dim = 3;
  std::vector<double> grad(dim * F.size());
  for (unsigned int ii = 0; ii<F.size(); ii++){
    for (int jj = 0; jj<dim; jj++){
      grad[dim*ii + jj] = -F[ii][jj];
    }
  }
  return grad;
}

void getEm_x(ElementMesh * em, double *x)
{
  int dim = em->dim;

  for (unsigned int ii = 0; ii<em->x.size(); ii++){
    for (int jj = 0; jj<dim; jj++){
      x[ii*dim + jj] = em->x[ii][jj] - em->X[ii][jj];
    }
  }
}

void setEm_x(ElementMesh * em, const double *x)
{
  int dim = em->dim;

  for (unsigned int ii = 0; ii<em->x.size(); ii++){
    for (int jj = 0; jj<dim; jj++){
      em->x[ii][jj] = em->X[ii][jj] + x[ii*dim + jj];
    }
  }
}

//world contains ElementMesh and LinSolver
//world also saves states to mesh files.
void solveStatics(World * world)
{
  ElementMesh * em = world->em_[0];
  int dim = em->dim;
  float E0 = potentialEnergy(em);
  LinSolver<double> * solver = world->solvers[0];
  std::vector<double> grad = staticGradient(em);
  grad.resize(grad.size() + 6,0);
  int iters = 0;
  double residual;
  
  std::vector<double> dx,x,x0;
  
  x0.resize(dim*em->x.size());
  dx.resize(x0.size() + 6);
  int ret = 0;
  int NEWTON_ITER = 5;
  double E1 = E0;
  for (int i = 0; i < NEWTON_ITER; i++){
    
    Eigen::SparseMatrix<double> K = em->getStiffnessSparse();
    fixTranslation(K, em);
    fixRotation(K, em);
    solver->solve(K, grad, dx, residual, iters);
    double h = 1;

    getEm_x(em, x0.data());
    ret = 0;
    //line search
    while (1){
      x = x0;
      //minimization steps in the negative gradient direction
      add_scaled(x, -h, dx);
      setEm_x(em, x.data());
      fem_error = 0;
      E1 = potentialEnergy(em);
      if ((E1>E0) || fem_error!= 0){
        h = 0.5f * h;
        if (h<1e-6){
          //std::cout << "SampleForce.cpp cannot find step size that reduces gradient norm\n";
          if (fem_error != 0){
            ret = -2;
          }
          break;
        }
      }
      else{
        break;
      }
    }
    if (ret < 0){
      break;
    }
    
    //world->saveMeshState();
    if (std::abs(E0 - E1) < 1e-4){
      break;
    }
    E0 = E1;

  }
}

//prints vertex positions and 
/// @param N grid resolution
void printPosEne(ElementMesh * em, std::ostream & out, int N)
{
  const int nV = 8;
  int vidx[nV];
  for (int i = 0; i < 2; i++){
    for (int j = 0; j < 2; j++){
      for (int k = 0; k < 2; k++){
        int linIdx = i * 4 + j * 2 + k;
        vidx[linIdx] = i*(N + 1)*(N + 1)*N + j*(N + 1)*N + k*N;
      }
    }
  }
  
  out << em->getEnergyElastic();
  for (int i = 0; i < nV; i++){
    out << " " << em->x[vidx[i]][0] << " " << em->x[vidx[i]][1] << " " << em->x[vidx[i]][2];
  }
  out << "\n";
}


void SampleForces(const ConfigFile & conf)
{
  int N = 4;
  World * world = new World();
  loadSamplingWorld(conf, world);
  ElementMesh * em = world->em_[0];
  FileUtilOut out(conf.dir+"/cube4.txt");
  world->em_[0]->saveMesh(out.out);
  out.close();

  world->simDirName = conf.dir + "/" + conf.getString("resultdir");

  std::string forceFile = conf.dir + "/" + conf.getString("forceList");
  std::string stretchForceFile = conf.dir + "/sforce.txt";
  std::vector<std::vector<float> > stretchForce = loadArr2d<float>(stretchForceFile);
  float forceScale = 1.0f;
  conf.getFloat("forceScale", forceScale);
  ForceSample forceSample;
  forceSample.forceScale = forceScale;

  forceSample.baseForce.resize(stretchForce.size());
  for (size_t i = 0; i < stretchForce.size(); i++){
    int nVert = 8;
    forceSample.baseForce[i].resize(nVert);
    for (int j = 0; j < nVert; j++){
      for (int k = 0; k < 3; k++){
        forceSample.baseForce[i][j][k] = forceScale * stretchForce[i][3 * j + k];
      }
    }
  }
  forceSample.ff = forceSample.baseForce;
  //forceSample.generateBase(); 
  //forceSample.scaleForces();
  
  int forceIdx = 0;
  conf.getInt("forceIdx", forceIdx);

  std::string combListFile = conf.dir + conf.getString("combList");
  std::vector<std::vector<int> > combList = loadArr2d<int>(combListFile);
  em->me.resize(em->e.size());

  for (size_t cidx = 0; cidx < 1; cidx++){
    for (size_t i = 0; i < em->e.size(); i++){
      em->me[i] = combList[cidx][i];
    }
    std::string dataFile = conf.dir +"/"+ std::to_string(cidx) + ".txt";
    FileUtilOut dataout(dataFile);
    dataout.out << forceSample.ff.size()<<" 25\n";
    for (size_t i = 0; i < forceSample.ff.size(); i++){
      copyForce(em, forceSample.ff[i]);
      //if (i% (2*forceSample.nMag) == 0){
      //  em->x = em->X;
      //}
      solveStatics(world);
      printPosEne(em, dataout.out, N);
      //if (i%forceSample.nMag == forceSample.nMag - 1){
      world->saveMeshState();
      world->frameCnt++;
      //}

      //  Eigen::Vector3d mn, mx;
      //  BBox(em->x, mn, mx);
      //  std::cout << "final bounding box\n";
      //  std::cout << mn[0] << " " << mn[1] << " " << mn[2] << "\n";
      //  std::cout << mx[0] << " " << mx[1] << " " << mx[2] << "\n";
  
      //float total = 0;
      //for (size_t fi = 0; fi < em->fe.size(); fi++){
      //  total += em->fe[fi].norm();
      //}
      //std::cout << "total force " << total << "\n";
    }
    dataout.close();
  }
}

///\@param p0 target for L-2 regularization
///\@param lambda regularization coeff
///\@param x output parameters
void leastSquares(const std::vector<std::vector<double > > & A,
  const std::vector<double> & b, const std::vector<double> & p0,
  double lambda,
  std::vector<double> & x)
{

}

void fitParams(const ConfigFile & conf)
{
  World * world = new World;
  world->simDirName = conf.getString("resultdir");
  std::string filename = conf.dir + "/" + conf.getString("coarseMeshFile");
  ElementMesh* em = new ElementMesh();
  world->em_.push_back(em);
  FileUtilIn in(filename);
  em->load(in.in);
  in.close();
  Quadrature q = Quadrature::NinePt;
  StrainEneANeo*ene = new StrainEneANeo();
  std::vector<Material *> material(1);
  MaterialQuad * matq = new MaterialQuad(ene, &q);
  for (unsigned int jj = 0; jj < matq->e.size(); jj++){
    matq->e[jj] = ene;
  }
  material[0] = matq;
  material[0]->density = 1e3;
  em->addMaterial(matq);
  for (size_t i = 0; i < em->me.size(); i++){
    em->me[i] = 0;
  }
  em->check();
  em->initElements(&q);
  int nV = em->e[0]->nV();
  int dim = 3;
  int nData = 0;

  Element * ele = em->e[0];
  //first two parameters are controlled by 1 parameter.
  int nParam = ene->param.size() - 1;
  std::vector<std::vector<double> > basisParam(nParam);
  basisParam[0] = std::vector<double>(ene->param.size(), 0);
  //E=1 nu=0.4
  basisParam[0][0] = 0.357143;
  basisParam[0][1] = 1.42857;
  for (int i = 1; i < nParam; i++){
    basisParam[i] = std::vector<double>(ene->param.size(), 0);
    basisParam[i][i + 1] = 1.0f;
  }
  world->frameCnt = 0;

  std::string dataDir = conf.dir + "/" + conf.getString("dataDir") + "/";
  int nComb = 1;
  conf.getInt("nComb", nComb);
  for (int combIdx = 0; combIdx < nComb; combIdx++){
    std::string datafile = dataDir + std::to_string(combIdx) + ".txt";
    FileUtilIn dataIn(datafile);
    dataIn.in >> nData;
    int nCol = 0;
    dataIn.in >> nCol;
    std::vector<double> Evec(nData, 0);
    std::vector<std::vector<Eigen::Vector3d> >  x(nData);
    for (int i = 0; i < nData; i++){
      dataIn.in >> Evec[i];
      x[i].resize(nV);
      for (int j = 0; j < nV; j++){
        for (int k = 0; k < 3; k++){
          dataIn.in >> x[i][j][k];
        }
      }
    }
    dataIn.close();
    std::vector<std::vector<double> > basisVal(nData);

    for (int i = 0; i < nData; i++){
      em->x = x[i];
      basisVal[i].resize(q.x.size() * basisParam.size());
      //world->saveMeshState();
      //world->frameCnt++;
      for (int j = 0; j < q.x.size(); j++){
        Eigen::Matrix3d F = ele->defGradCached(j, em->X, em->x);
        for (size_t k = 0; k < basisParam.size(); k++){
          matq->e[j]->param = basisParam[k];
          double E = ele->detJ[j] * q.w[j] * matq->e[j]->getEnergy(F);
          basisVal[i][j*basisParam.size() + k] = E;
        }
      }
    }

    std::string basisOutFile = conf.dir + "/basis" + std::to_string(combIdx) + ".txt";
    FileUtilOut basisOut(basisOutFile);
    basisOut.out << nData << " " << basisVal[0].size() << "\n";
    for (int i = 0; i < nData; i++){
      for (size_t j = 0; j < basisVal[i].size(); j++){
        basisOut.out << basisVal[i][j] << " ";
      }
      basisOut.out << "\n";
    }
    basisOut.close();
  }

  //double lambda = 1.0f;
  //std::vector<double> p;
  //std::vector<double> p0;

  //std::string combListFile = conf.dir + conf.getString("combList");
  //std::vector<std::vector<int> > combList = loadArr2d<int>(combListFile);
  //int combIdx = 0;
  //conf.getInt("combIdx", combIdx);

  //leastSquares(basisVal, Evec, p0, lambda, p);
}

std::vector<int> extractOneComb(int i0, int j0, int k0,
  const std::vector<int> & gridres, 
  const std::vector<int> & a,
  int coarseSize)
{
  int cnt = 0;
  std::vector <int> comb(coarseSize*coarseSize*coarseSize,0);
  for (int i = i0; i < i0 + coarseSize; i++){
    for (int j = j0; j < j0 + coarseSize; j++){
      for (int k = k0; k < k0 + coarseSize; k++){
        int idx = gridToLinearIdx(i, j, k, gridres);
        int mat = a[idx];
        comb[cnt] = mat;
        cnt++;
      }
    }
  }
  return comb;
}

void assignOneComb(int i0, int j0, int k0,
  const std::vector<int> & gridres,
  std::vector<int> & a,
  int coarseSize,
  const std::vector<int> comb)
{
  int cnt = 0;
  for (int i = i0; i < i0 + coarseSize; i++){
    for (int j = j0; j < j0 + coarseSize; j++){
      for (int k = k0; k < k0 + coarseSize; k++){
        int mat = comb[cnt];
        int idx = gridToLinearIdx(i, j, k, gridres);
        a[idx] = mat;
        cnt++;
      }
    }
  }
}

void extractComb(const ConfigFile & conf)
{
  int dim = 3;
  conf.getInt("dim", dim);
  ElementMesh * em = new ElementMesh();
  em->dim = dim;
  std::string meshFile = conf.getString("meshfile");
  FileUtilIn in(conf.dir + "/" + meshFile);
  em->load(in.in);

  //size of a coarse element block
  int coarseSize = 4;
  conf.getInt("coarseSize", coarseSize);
  std::vector<int> gridres;
  conf.getIntVector("gridres", gridres);
  std::map < std::vector<int>, int> matCombs;
  int matID = 0;
  
  ///resolution of the coarse mesh
  std::vector<int> coarseRes(3, 0);
  for (size_t i = 0; i < gridres.size(); i++){
    coarseRes[i] = gridres[i] / coarseSize;
  }
  std::vector<int> coarseMat(coarseRes[0] * coarseRes[1] * coarseRes[2],0);

  for (int i = 0; i <= gridres[0]-coarseSize; i+=coarseSize){
    for (int j = 0; j <= gridres[1] - coarseSize; j += coarseSize){
      for (int k = 0; k <= gridres[2] - coarseSize; k += coarseSize){
        std::vector<int> comb = extractOneComb(i, j, k, gridres, em->me, coarseSize);
        auto it = matCombs.find(comb);
        int coarseIdx = gridToLinearIdx(i / coarseSize, j / coarseSize, k / coarseSize, coarseRes);
        if (it == matCombs.end()){
          matCombs[comb] = matID;
          coarseMat[coarseIdx] = matID;
          matID++;          
        }
        else{
          coarseMat[coarseIdx] = it->second;
        }
      }
    }
  }
  std::cout << matID << " combinations.\n";
  std::string combFile = conf.dir + "/" + conf.getString("combFile");
  FileUtilOut out(combFile);
  int nCol = coarseSize*coarseSize*coarseSize;
  out.out << matID << " " << nCol << "\n";
  std::vector<std::vector< int > > combList(matID);
  for (auto it = matCombs.begin(); it != matCombs.end(); it++){
    int id = it->second;
    combList[id] = it->first;
  }
  for (size_t i = 0; i < combList.size(); i++){
    for (size_t j = 0; j < combList[i].size(); j++){
      out.out << combList[i][j] << " ";
    }
    out.out << "\n";
  }
  out.close();

  //save coarse mesh.
  ElementRegGrid * coarseMesh = new ElementRegGrid(coarseRes[0], coarseRes[1], coarseRes[2]);
  for (size_t i = 0; i < coarseMesh->me.size(); i++){
    coarseMesh->me[i] = coarseMat[i];
  }
  std::string coarseFile = conf.dir + "/coarseFEM.txt";
  FileUtilOut coarseOut(coarseFile);
  coarseMesh->saveMesh(coarseOut.out);
  coarseOut.close();

  //reconstruct rigid material for debugging and validation.
  bool debugging = false;
  if (debugging){
    for (size_t i = 0; i < em->me.size(); i++){
      em->me[i] = 0;
    }
    for (int i = 0; i <= gridres[0] - coarseSize; i += coarseSize){
      for (int j = 0; j <= gridres[1] - coarseSize; j += coarseSize){
        for (int k = 0; k <= gridres[2] - coarseSize; k += coarseSize){
          int coarseIdx = gridToLinearIdx(i / coarseSize, j / coarseSize, k / coarseSize, coarseRes);
          int coarseMatId = coarseMesh->me[coarseIdx];
          std::vector<int> comb = combList[coarseMatId];
          assignOneComb(i, j, k, gridres, em->me, coarseSize, comb);
        }
      }
    }
    savePartObj(em, 0, conf.dir + "/debugMat0.obj");
  }
  
}

//load coarse materials
//with 5 quadrature points and
//StrainEneANeo
void makeCoarseMats(std::vector<std::vector<double > > & params, 
  std::vector<Material*>& materials)
{
  //hard coded poisson's ratio. 
  //load from file later.
  //double nu = 0.4;
  double density = 1e3;
  double basisParam[2] = { 0.357143, 1.42857 };
  //#params per quadrature point.
  int nParam = params[0].size();
  
  Quadrature * q = new Quadrature();
  *q = Quadrature::NinePt;
  materials.resize(params.size());
  for (size_t i = 0; i < params.size(); i++){
    MaterialQuad * mat = new MaterialQuad();
    int nParamPerQuad = nParam / mat->e.size();
    mat->q = q;
    mat->e.resize(q->x.size());
    mat->density = density;
    for (size_t j = 0; j < mat->e.size(); j++){
      StrainEneANeo * se = new StrainEneANeo();
      se->param[0] = basisParam[0] * params[i][nParamPerQuad * j];
      se->param[1] = basisParam[1] * params[i][nParamPerQuad * j];
      for (int k = 1; k < nParamPerQuad; k++){
        se->param[k + 1] = params[i][nParamPerQuad * j + k];
      }
      mat->e[j] = se;
    }
    materials[i] = mat;
  }
}
