#include "StretchMesh.hpp"

#include "ArrayUtil.hpp"
#include "Element.hpp"

#include <iomanip>

void ind2sub(int idx, int N, int & i, int & j)
{
  i = idx / N;
  j = idx - i*N;
}

void ind2sub(int idx, int N1, int N2, int & i, int & j, int &k)
{
  i = idx / (N1*N2);
  idx -= i*N1*N2;
  j = idx/N2;
  k = idx - j*N2;
}

int sub2ind(int i, int j, int N)
{
  return i*N+j;
}

int sub2ind(int i, int j, int k, int N1, int N2)
{
  return i*N1*N2+j*N2+k;
}

void stretchXRight(ElementMesh & em, const StretchOpts & opts,
                 int targetLabel)
{
  float gridlen = (float)(em.X[em.e[0]->at(1)][1] - em.X[em.e[0]->at(0)][1]);
  double offset = opts.offset[2];
  float linearScale = (float)(offset / opts.nMidEle);
  std::vector<int> stretched(em.X.size(), 0);
  std::vector<Eigen::Vector3d> X0 = em.X;
  float m = opts.xmax + opts.beamwidth + gridlen;
  for(unsigned int ii = 0; ii<em.e.size(); ii++){
    Eigen::Vector3d v0 = X0[em.e[ii]->at(0)];
    Eigen::Vector3d v2 = X0[em.e[ii]->at(2)];
    if(opts.label[ii] == targetLabel && v2[0]<m){
      continue;
    }else if(opts.label[ii] == targetLabel && v0[0]>m){
      int steps = (int)((v2[0]-m) / gridlen);
      float dx_right = (steps+1) * linearScale;
      float dx_left =  (steps)    * linearScale;
      for(int jj = 0; jj<em.e[ii]->nV(); jj++){
        int vidx = em.e[ii]->at(jj);
        if(stretched[vidx]){
          continue;
        }
        if(jj<2){
          em.X[vidx][0] += dx_left;
        }else{
          em.X[vidx][0] += dx_right;
        }
        stretched[vidx] = 1;
      }
    }
  }
}

void stretchX(ElementMesh & em, const StretchOpts & opts,
                 int targetLabel)
{
  float gridlen = (float)(em.X[em.e[0]->at(1)][1] - em.X[em.e[0]->at(0)][1]);
  double offset = -opts.offset[0];
//  if(targetLabel==2){
//    offset = -opts.offset[2];
//  }
  double linearScale = offset/opts.nMidEle;
  std::vector<int> stretched(em.X.size(), 0);
  std::vector<Eigen::Vector3d> X0 = em.X;
  for(unsigned int ii = 0; ii<em.e.size(); ii++){
    Eigen::Vector3d v0 = X0[em.e[ii]->at(0)];
    Eigen::Vector3d v2 = X0[em.e[ii]->at(2)];
    if(opts.label[ii] == targetLabel && v2[0]<opts.xmin){
      for(int jj = 0; jj<em.e[ii]->nV(); jj++){
        int vidx = em.e[ii]->at(jj);
        if(stretched[vidx]){
          continue;
        }
        em.X[vidx][0] += offset;
        stretched[vidx] = 1;
      }
    }else if(opts.label[ii] == targetLabel && v2[0]<opts.xmax){
      int steps = (int)((opts.xmax-v0[0]) / gridlen);
      float dx_right = (steps) * linearScale;
      float dx_left =  (steps+1)    * linearScale;
      for(int jj = 0; jj<em.e[ii]->nV(); jj++){
        int vidx = em.e[ii]->at(jj);
        if(stretched[vidx]){
          continue;
        }
        if(jj<2){
          em.X[vidx][0] += dx_left;
        }else{
          em.X[vidx][0] += dx_right;
        }
        stretched[vidx] = 1;
      }
    }
  }

}

void stretchY(ElementMesh & em, const StretchOpts & opts)
{
  float gridlen = em.X[em.e[0]->at(1)][1] - em.X[em.e[0]->at(0)][1];
  double offset = opts.offset[1];
  double linearScale = offset/opts.nBeamEle;
  std::vector<int> stretched(em.X.size(), 0);
  for(unsigned int ii = 0; ii<em.e.size(); ii++){
    int label = opts.label[ii];
    if(label==0){
      for(int jj = 0; jj<em.e[ii]->nV(); jj++){
        int vidx = em.e[ii]->at(jj);
        if(stretched[vidx]){
          continue;
        }
        em.X[vidx][1] += offset;
        stretched[vidx] = 1;
      }
    }else if(label == 1){
      for(int jj = 0; jj<em.e[ii]->nV(); jj++){
        int vidx = em.e[ii]->at(jj);
        if(stretched[vidx]){
          continue;
        }
        double dx = linearScale * (em.X[vidx][1]-opts.beamy0)/gridlen;
//        std::cout<<ii<<" "<<dx<<" "<<em.X[vidx][1]<<"\n";
        em.X[vidx][1] += dx;
        stretched[vidx] = 1;
      }
    }
  }
}

void
stretchJumper(ElementMesh & em, const StretchOpts & opts)
{
//  stretchX(em, opts,0);
//  stretchX(em, opts,2);
  stretchY(em, opts);
  stretchXRight(em, opts, 2);
  em.x = em.X;
}

void printGrid(const std::vector<std::vector< int> > & grid,
               const std::vector<int> & label)
{
  int nx = (int)grid.size();
  int ny = (int)grid[0].size();
  for(int jj = ny-1; jj>=0; jj--){
    for(int ii = 0; ii<nx; ii++){
      std::cout<<std::setfill(' ') << std::setw(3);
      if(grid[ii][jj]>=0){
        std::cout<<(1+label[grid[ii][jj]])<< " ";
      }else{
        std::cout<<0<< " ";
      }
    }
    std::cout<<"\n";
  }
}

void printGrid3d(std::vector<std::vector< std::vector<int> > > & grid)
{
  int nx = (int)grid.size();
  int ny = (int)grid[0].size();
  int nz = (int)grid[0][0].size();
  for(int kk = 0; kk<nz; kk++){
  for(int jj = ny-1; jj>=0; jj--){
    for(int ii = 0; ii<nx; ii++){
      std::cout<<std::setfill(' ') << std::setw(3);
      std::cout<<grid[ii][jj][kk]<< " ";
    }
    std::cout<<"\n";
  }
  std::cout<<"\n";
  }
}

void
segmentJumper(const ElementMesh & em, StretchOpts & opts)
{
  //assume input mesh em has uniform grid.
  float gridlen = em.X[em.e[0]->at(1)][1] - em.X[em.e[0]->at(0)][1];
  std::vector<std::vector< int> > grid;
  em.grid2D(grid);
  int nx = (int)grid.size();
  int ny = (int)grid[0].size();
  //find top right element
  int topj = -1;
  int topi = -1;
  for(int jj = ny-1 ; jj>=0; jj-- ){
    for( int ii = nx-1; ii>=0; ii--){
      if(grid[ii][jj]>=0){
        topi = ii;
        topj = jj;
        break;
      }
    }
    if(topi>=0){
      break;
    }
  }

  //find top left vertex
  int ii = topi-1;
  for(; ii>=0; ii--){
    if(grid[ii][topj]<0){
      break;
    }
  }
  ii++;
  ii+=(int)(opts.leftmargin/gridlen+0.5);
  int xminIdx = ii, xmaxIdx = -1;

  int beamstart = -1;
  int beamend = -1;

  float thresh=0.5*gridlen + opts.beamwidth;
  for(int jj = topj; jj>=0; jj--){
    if(grid[topi][jj]<0){
      //probably unexpected model.
      break;
    }
    float beamWidth = gridlen;
    int ii = topi-1;
    for(; ii>=0; ii--){
      if(grid[ii][jj]>=0){
        beamWidth += gridlen;
      }else{
        break;
      }
    }
    if(beamstart<0 && beamWidth < thresh){
      xmaxIdx = ii;
      beamstart = jj;
    }
    if(beamstart >= 0 && beamWidth>thresh){
      beamend = jj+1;
      break;
    }
  }

  if(beamstart < 0 || beamend<0){
    std::cout<<"Error segment fem. Cannot find beam.";
  }
  int tmp = beamend;
  beamend = beamstart;
  beamstart = tmp;
//  std::cout<<"beam start "<<beamstart<<" "<<beamend<<"\n";
  opts.label.resize(em.e.size(), -1);

  //label beam as 1
  for(int jj = beamstart; jj<=beamend; jj++){
    for(int ii=topi; ii>=0; ii--){
      int eIdx = grid[ii][jj];
      if(eIdx<0){
        break;
      }
      opts.label[eIdx] = 1;
    }
  }

  //flood top
  std::vector<int> q;
  int idx = sub2ind(topi, topj, ny);
  q.push_back(idx);
  const int NNBR=4;
  int NBRS[NNBR][2] = {{-1,0},{1,0},{0,-1},{0,1}};
  while(!q.empty()){
    idx = q.back();
    q.pop_back();
    int ii, jj;
    ind2sub(idx, ny, ii, jj);
    int eIdx = grid[ii][jj];
    opts.label[eIdx] = 0;
    for(int nbr = 0; nbr<NNBR; nbr++){
      int ni = ii + NBRS[nbr][0];
      int nj = jj + NBRS[nbr][1];
      if(ni<0 || ni>=nx || nj<0 || nj>=ny){
        continue;
      }
      //neighbor element
      int ne=grid[ni][nj];
      if(ne<0){
        continue;
      }
      if(opts.label[ne]>=0){
        continue;
      }
      q.push_back(sub2ind(ni,nj,ny));
    }
  }

  for(unsigned int ii = 0; ii<opts.label.size(); ii++){
    if(opts.label[ii]<0){
      opts.label[ii] = 2;
    }
  }
  //printGrid(grid,opts.label);

//  std::cout<<"horizontal range "<<xminIdx<<" "<<xmaxIdx<<"\n";

  int eidx = grid[xminIdx][topj];
  opts.xmin = 0.5 * (em.X[em.e[eidx]->at(2)][0] + em.X[em.e[eidx]->at(0)][0]);
  eidx = grid[xmaxIdx][topj];
  opts.xmax = 0.5 * (em.X[em.e[eidx]->at(2)][0] + em.X[em.e[eidx]->at(0)][0]);

  opts.nMidEle = xmaxIdx - xminIdx + 1;
//  std::cout<<"number of stretched elements " << opts.nMidEle<<"\n";
  opts.nBeamEle = beamend - beamstart+1;
  eidx = grid[topi][beamstart];
  opts.beamy0 = em.X[em.e[eidx]->at(0)][1];
//  std::cout<<"beam y0 "<<opts.beamy0<<"\n";
}

void check_square_ele(const ElementMesh * em)
{
  if(em->dim != 3){
    return;
  }
  double eps = 1e-6;
  for(int ei = 0; ei<em->e.size(); ei++){
    Element* ele = em->e[ei];
    Eigen::Vector3d eSize = em->X[ele->at(7)] - em->X[ele->at(0)];
    for(int ii = 1; ii<ele->nV()-1; ii++){
      Eigen::Vector3d disp = em->X[ele->at(ii)] - em->X[ele->at(0)];
      for(int jj = 0; jj<em->dim; jj++){
        double val = std::abs(disp[jj]);
        if(! (val<eps || std::abs(val - eSize[jj])<eps)){
          std::cout<<ei<<" "<<ii<<"\n";
        }
      }
    }
  }
}

void
stretchJumper3d(ElementMesh & em, const StretchOpts & opts)
{
  //assume input mesh em has uniform grid.
  float gridlen = em.X[em.e[0]->at(4)][0] - em.X[em.e[0]->at(0)][0];
  //First, stretch horizontal.
  //apply linear scaling to middle vertices
  double linearScale = opts.offset[0]/opts.nMidEle;
//  std::cout<<"Stretch "<<opts.offset[0]<<" "<<opts.offset[1]<<"\n";
//  std::cout<<"linear Scale "<< linearScale <<"\n";
//  std::cout<<opts.nMidEle<<" "<<opts.xmin<<" "<<opts.xmax<<"\n";
  std::vector<int> stretched(em.X.size(), 0);
  std::vector<Eigen::Vector3d> X0 = em.X;
  //stretch toward right
  for(unsigned int ii = 0; ii<em.e.size(); ii++){
    Eigen::Vector3d v0 = X0[em.e[ii]->at(0)];
    Eigen::Vector3d v2 = X0[em.e[ii]->at(2)];
    if(v2[0]<opts.xmin){
      continue;
    }else if(v0[0]<opts.xmax){
      int steps = (int)((v2[0] - opts.xmin) / gridlen);
      float dx_right = (steps+1) * linearScale;
      float dx_left =   steps    * linearScale;
      for(int jj = 0; jj<em.e[ii]->nV(); jj++){
        int vidx = em.e[ii]->at(jj);
        if(stretched[vidx]){
          continue;
        }
        if(jj<4){
          em.X[vidx][0] += dx_left;
        }else{
          em.X[vidx][0] += dx_right;
        }
        stretched[vidx] = 1;
      }
    }else{
      for(int jj = 0; jj<em.e[ii]->nV(); jj++){
        int vidx = em.e[ii]->at(jj);
        if(stretched[vidx]){
          continue;
        }
        em.X[vidx][0] += opts.offset[0];
        stretched[vidx] = 1;
      }
    }
  }

  //stretch vertical upwards
  linearScale = opts.offset[1]/opts.nBeamEle;
  std::fill(stretched.begin(), stretched.end(), 0);
  for(unsigned int ii = 0; ii<em.e.size(); ii++){
    int label = opts.label[ii];
    if(label==0){
      for(int jj = 0; jj<em.e[ii]->nV(); jj++){
        int vidx = em.e[ii]->at(jj);
        if(stretched[vidx]){
          continue;
        }
        em.X[vidx][1] += opts.offset[1];
        stretched[vidx] = 1;
      }
    }else if(label == 1){
      for(int jj = 0; jj<em.e[ii]->nV(); jj++){
        int vidx = em.e[ii]->at(jj);
        if(stretched[vidx]){
          continue;
        }
        double dx = linearScale * (em.X[vidx][1]-opts.beamy0)/gridlen;
//        std::cout<<ii<<" "<<dx<<" "<<em.X[vidx][1]<<"\n";
        em.X[vidx][1] += dx;
        stretched[vidx] = 1;
      }
    }
  }
  em.x = em.X;

//  check_square_ele(&em);
}

void
segmentJumper3d(const ElementMesh & em, StretchOpts & opts)
{
  //assume input mesh em has uniform grid.
  float gridlen = em.X[em.e[0]->at(4)][0] - em.X[em.e[0]->at(0)][0];
  std::vector<std::vector< std::vector<int> > > grid;
  em.grid3D(grid);
  int nx = (int)grid.size();
  int ny = (int)grid[0].size();
  int nz = (int)grid[0][0].size();
  //find top right element
  int topj = -1;
  int topi = -1;
  int topk = -1;
  for (int ii = nx - 1; ii >= 0; ii--){
    for(int jj = ny-1 ; jj>=ny/2; jj-- ){
      for(int kk = nz-1; kk>=0; kk--){
        if(grid[ii][jj][kk]>=0){
          topi = ii;
          topj = jj;
          topk = kk;
          break;
        }
      }
      if(topk>=0){
        break;
      }
    }
    if(topi>=0){
      break;
    }
  }

  //find top left vertex
  int ii;
  for (ii = topi - 1; ii >= 0; ii--){
    if(grid[ii][topj][topk]<0){
      break;
    }
  }
  ii++;
  ii+=(int)(opts.leftmargin/gridlen+0.5);
  int xminIdx = ii, xmaxIdx = -1;

//  printGrid3d(grid);
//  std::cout<<"Top ele "<< topi<<" "<<topj<<"\n";

  int beamstart = -1;
  int beamend = -1;

  float thresh=0.5*gridlen + opts.beamwidth;
  for(int jj = topj; jj>=0; jj--){
    if(grid[topi][jj][topk]<0){
      //probably unexpected model.
      break;
    }
    float beamWidth = gridlen;
    int ii;
    for (ii = topi - 1; ii >= 0; ii--){
      if(grid[ii][jj][topk]>=0){
        beamWidth += gridlen;
      }else{
        break;
      }
      if(beamWidth>thresh){
        break;
      }
    }
    if(beamstart<0 && beamWidth < thresh){
      xmaxIdx = ii;
      beamstart = jj;
    }
    if(beamstart >= 0 && (beamWidth>thresh)){
      beamend = jj+1;
      break;
    }
    if (beamstart >= 0 && (beamWidth <= thresh && jj == 0)) {
      beamend = jj;
      break;
    }
  }

  if(beamstart < 0 || beamend<0){
    std::cout<<"Error segment fem. Cannot find beam.";
  }
  int tmp = beamend;
  beamend = beamstart;
  beamstart = tmp;
//  std::cout<<"beam start "<<beamstart<<" "<<beamend<<"\n";
  opts.label.resize(em.e.size(), -1);

  //label beam as 1
  for(int jj = beamstart; jj<=beamend; jj++){
    for(int ii=topi; ii>=0; ii--){
      int eIdx = grid[ii][jj][topk];
      if(eIdx<0){
        break;
      }
      for(int kk = topk; kk>=0; kk--){
        int eIdx = grid[ii][jj][kk];
        if(eIdx<0){
          continue;
        }
        opts.label[eIdx] = 1;
      }
    }
  }

  //flood top
  std::vector<int> q;
  int idx = sub2ind(topi, topj, topk, ny, nz);
  q.push_back(idx);
  const int NNBR=6;
  int NBRS[NNBR][3] = {{-1,0,0},{1,0,0},{0,-1,0},{0,1,0},{0,0,-1},{0,0,1}};
  while(!q.empty()){
    idx = q.back();
    q.pop_back();
    int ii, jj,kk;
    ind2sub(idx, ny,nz, ii, jj,kk);
    int eIdx = grid[ii][jj][kk];
    opts.label[eIdx] = 0;
    for(int nbr = 0; nbr<NNBR; nbr++){
      int ni = ii + NBRS[nbr][0];
      int nj = jj + NBRS[nbr][1];
      int nk = kk + NBRS[nbr][2];
      if(ni<0 || ni>=nx || nj<0 || nj>=ny || nk<0 || nk>=nz){
        continue;
      }
      //neighbor element
      int ne=grid[ni][nj][nk];
      if(ne<0){
        continue;
      }
      if(opts.label[ne]>=0){
        continue;
      }
      q.push_back(sub2ind(ni,nj,nk,ny,nz));
    }
  }

  for(unsigned int ii = 0; ii<opts.label.size(); ii++){
    if(opts.label[ii]<0){
      opts.label[ii] = 2;
    }
  }

  //print segmentation to console
//  std::vector<std::vector<std::vector<int> > > labels = grid;
//  for(int ii = 0; ii<nx; ii++){
//    for(int jj = 0; jj<ny; jj++){
//      for(int kk = 0; kk<nz; kk++){
//        int eidx = grid[ii][jj][kk];
//        if (eidx<0){
//          continue;
//        }
//        labels[ii][jj][kk] = opts.label[eidx];
//      }
//    }
//  }
//  printGrid3d(labels);

//  std::cout<<"horizontal range "<<xminIdx<<" "<<xmaxIdx<<"\n";

  int eidx = grid[xminIdx][topj][topk];
  opts.xmin = 0.5 * (em.X[em.e[eidx]->at(2)][0] + em.X[em.e[eidx]->at(0)][0]);
  eidx = grid[xmaxIdx][topj][topk];
  opts.xmax = 0.5 * (em.X[em.e[eidx]->at(2)][0] + em.X[em.e[eidx]->at(0)][0]);

  opts.nMidEle = xmaxIdx - xminIdx + 1;
//  std::cout<<"number of stretched elements " << opts.nMidEle<<"\n";
  opts.nBeamEle = beamend - beamstart+1;
  eidx = grid[topi][beamstart][topk];
  opts.beamy0 = em.X[em.e[eidx]->at(0)][1];
//  std::cout<<"beam y0 "<<opts.beamy0<<"\n";
}
