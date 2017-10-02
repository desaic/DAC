/*************************************************************************\

  Copyright 2012 The University of North Carolina at Chapel Hill.
  All Rights Reserved.

  Permission to use, copy, modify and distribute this software and its
  documentation for educational, research and non-profit purposes, without
   fee, and without a written agreement is hereby granted, provided that the
  above copyright notice and the following three paragraphs appear in all
  copies.

  IN NO EVENT SHALL THE UNIVERSITY OF NORTH CAROLINA AT CHAPEL HILL BE
  LIABLE TO ANY PARTY FOR DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR
  CONSEQUENTIAL DAMAGES, INCLUDING LOST PROFITS, ARISING OUT OF THE
  USE OF THIS SOFTWARE AND ITS DOCUMENTATION, EVEN IF THE UNIVERSITY
  OF NORTH CAROLINA HAVE BEEN ADVISED OF THE POSSIBILITY OF SUCH
  DAMAGES.

  THE UNIVERSITY OF NORTH CAROLINA SPECIFICALLY DISCLAIM ANY
  WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
  MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.  THE SOFTWARE
  PROVIDED HEREUNDER IS ON AN "AS IS" BASIS, AND THE UNIVERSITY OF
  NORTH CAROLINA HAS NO OBLIGATIONS TO PROVIDE MAINTENANCE, SUPPORT,
  UPDATES, ENHANCEMENTS, OR MODIFICATIONS.

  The authors may be contacted via:

  US Mail:       GAMMA Research Group at UNC
                       Department of Computer Science
                       Sitterson Hall, CB #3175
                       University of N. Carolina
                       Chapel Hill, NC 27599-3175

  Phone:            (919)962-1749

  EMail:              geom@cs.unc.edu; tang_m@zju.edu.cn


\**************************************************************************/

// Author: Tang, Min tang_m@zju.edu.cn
#include <iostream>
#include "DeformModel.h"

#include "timing.h"
extern CBVHTimer tm;

#ifndef swapI
#define swapI(a, b) {\
	unsigned int tmp = a;\
	a = b;\
	b = tmp;\
}
#endif

extern float
Intersect_VF(const vec3f &ta0, const vec3f &tb0, const vec3f &tc0,
			 const vec3f &ta1, const vec3f &tb1, const vec3f &tc1,
			 const vec3f &q0, const vec3f &q1,
			 vec3f &qi, vec3f &baryc);
extern float
Intersect_EE(const vec3f &ta0, const vec3f &tb0, const vec3f &tc0, const vec3f &td0,
			 const vec3f &ta1, const vec3f &tb1, const vec3f &tc1, const vec3f &td1,
			 vec3f &qi);

float
DeformModel::intersect_vf(unsigned int fid1, unsigned int vid2, unsigned int fid2)
{
	if (!_fac_boxes[fid1].overlaps(_vtx_boxes[vid2]))
		return -1.f;

	for (id_list::iterator it1=_vtx_fids[vid2].begin(); it1!=_vtx_fids[vid2].end(); it1++) {
			unsigned int fid = *it1;
			if (!Covertex_F(fid1, fid)) {
				if (fid == fid2)
					return do_vf(fid1, vid2);
				else
					return -1.f;
			}
	}
	return -1.f;
}

float
DeformModel::intersect_vf(unsigned int fid, unsigned int vid)
{
	if (!_fac_boxes[fid].overlaps(_vtx_boxes[vid]))
		return -1.f;

	return do_vf(fid, vid);
}

float
DeformModel::do_vf(unsigned int fid, unsigned int vid)
{
	_counts[omp_get_thread_num()]._num_vf_test++;

	vec3f qi, baryc;
	unsigned v0 = _tris[fid].id0();
	unsigned v1 = _tris[fid].id1();
	unsigned v2 = _tris[fid].id2();

	//tm.startTiming(8);

	float ret = Intersect_VF(
		_prev_vtxs[v0],  _prev_vtxs[v1], _prev_vtxs[v2],
		_cur_vtxs[v0],   _cur_vtxs[v1],  _cur_vtxs[v2],
		_prev_vtxs[vid], _cur_vtxs[vid], qi, baryc);
	//tm.endTiming(8);

	if (ret> -0.5) {
    //std::cout << "vf " <<fid << " " << vid << "\n";
    Contactvf vf;
    vf.f = fid;
    vf.v = vid;
    vflist.push_back(vf);
		_counts[omp_get_thread_num()]._num_vf_true++;
	}

	return ret;
}

float
DeformModel::intersect_ee(unsigned int e1, unsigned int e2, unsigned int f1, unsigned int f2)
{
	if (!_edg_boxes[e1].overlaps(_edg_boxes[e2]))
		return -1.f;

	unsigned int e[2];
	unsigned int f[2];

	if (e1 > e2) {
		e[0] = e1, e[1] = e2;
		f[0] = f1, f[1] = f2;
	} else {
		e[0] = e2, e[1] = e1;
		f[0] = f2, f[1] = f1;
	}

	for (int i=0; i<2; i++)
		for (int j=0; j<2; j++) {
			unsigned int ff1 = _edges[e[0]].fid(i);
			unsigned int ff2 = _edges[e[1]].fid(j);

			if (ff1 == -1 || ff2 == -1)
				continue;

			if (!Covertex_F(ff1, ff2)) {
					if (ff1 == f[0] && ff2 == f[1])
						return do_ee(e1, e2);
					else
						return -1.f;
			}
		}

	return -1.f;
}

float
DeformModel::intersect_ee(unsigned int e1, unsigned int e2)
{
	if (!_edg_boxes[e1].overlaps(_edg_boxes[e2]))
		return -1.f;

	return do_ee(e1, e2);
}

float
DeformModel::do_ee(unsigned int e1, unsigned int e2)
{
	_counts[omp_get_thread_num()]._num_ee_test++;

	vec3f qi;
	unsigned v0 = _edges[e1].vid(0);
	unsigned v1 = _edges[e1].vid(1);
	unsigned w0 = _edges[e2].vid(0);
	unsigned w1 = _edges[e2].vid(1);


  //if (v0 == 23 && v1 == 10 && w0 == 18 && w1 == 17) {
  //  std::cout << "debug!\n";
  //}
	//tm.startTiming(7);
	float ret = Intersect_EE(
		_prev_vtxs[v0], _prev_vtxs[v1],
		_prev_vtxs[w0], _prev_vtxs[w1],
		_cur_vtxs[v0], _cur_vtxs[v1],
		_cur_vtxs[w0], _cur_vtxs[w1],
		qi);
	//tm.endTiming(7);

	if (ret> -0.5) {
    //std::cout << "ee " << v0 << " " << v1 << "\n";
    //std::cout << "   " << w0 << " " << w1 << "\n";
    Contactee ee;
    ee.e0[0] = v0;
    ee.e0[1] = v1;
    ee.e1[0] = w0;
    ee.e1[1] = w1;
    eelist.push_back(ee);
		_counts[omp_get_thread_num()]._num_ee_true++;
	}

	return ret;
}

void
DeformModel::test_feature_0(unsigned id1, unsigned int id2)
{
	// 6 VF test
	for (int i=0; i<3; i++) {
		intersect_vf(id1, _tris[id2].id(i), id2);
		intersect_vf(id2, _tris[id1].id(i), id1);
	}

	// 9 EE test
	for (int i=0; i<3; i++)
		for (int j=0; j<3; j++) {
		unsigned int e0 = _tri_edges[id1].id(i);
		unsigned int e1 = _tri_edges[id2].id(j);
		
		intersect_ee(e0, e1, id1, id2);
	}
}
