#ifndef STRETCH_MESH_HPP
#define STRETCH_MESH_HPP

#include "TrigMesh.hpp"
#include "ElementMesh.hpp"
#include <Eigen/Dense>
#include <vector>

struct StretchOpts
{
  Eigen::Vector3d offset;
  ///@brief 0 top. 1 middle beam. 2 bottom.
  std::vector<int> label;
  ///@brief interval for horizontal stretching
  float xmin, xmax;
  ///@brief maximum width of the bending beam
  float beamwidth;
  float leftmargin;
  ///@brief number of elements in the middle that can be scaled.
  int nMidEle;

  ///@brief length of beam measured in number of elements.
  int nBeamEle;
  ///@brief lowest point of the beam.
  float beamy0;
  StretchOpts():offset(Eigen::Vector3d::Zero()),
  xmin(0),xmax(0),
  beamwidth(0),
  leftmargin(0),
  nMidEle(0),
  nBeamEle(0),
  beamy0(0){}
};

///@brief label and range is written into opts
/// Top = 0
/// Beam = 1
/// Bottom = 2
void segmentJumper(const ElementMesh & em, StretchOpts & opts);

///@brief stretch jumper according to labels and options in opts.
void stretchJumper(ElementMesh &em, const StretchOpts & opts);

void
stretchJumper3d(ElementMesh & em, const StretchOpts & opts);

void
segmentJumper3d(const ElementMesh & em, StretchOpts & opts);

#endif
