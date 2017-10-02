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

#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#include <gl/gl.h>

#include <math.h>
#include <stdio.h>
#include <assert.h>

#include "timing.h"
extern CBVHTimer tm;

#include "tri_pair.h"
#include <vector>
#include <algorithm>
using namespace std;

non_adjacent_pair_list non_adj_list;

#include "bvh_front.h"
bvh_front_list front_list;

#include "DeformModel.h"

#include "UpperBVH.h"

#include <omp.h>

DeformModel::DeformModel()
{
	Init();
}

DeformModel::~DeformModel()
{
	Clean();
}

void
DeformModel::Init()
{
	_num_vtx = 0;
	_num_frame = 0;

	_vtxs = NULL;
	_cur_vtxs = NULL;
	_prev_vtxs = NULL;
	_vtx_fids = NULL;

	_colors = NULL;
	_nrms = NULL;

	_num_tri = 0;
	_tris = NULL;
	_tri_edges = NULL;
	_tri_nrms = NULL;
	_old_tri_nrms = NULL;
	_tri_flags = NULL;

	_num_edge = 0;
	_edges = NULL;

	_tree = NULL;
	_tri_centers = NULL;
	_tri_boxes = NULL;
	
	for (int i=0; i<16; i++)
		_counts[i].init();

	_num_parts = 0;
	_parts = NULL;
}

void
DeformModel::Clean()
{
	if (_vtxs) delete [] _vtxs;
	if (_vtx_fids) delete [] _vtx_fids;
	if (_cur_vtxs) delete [] _cur_vtxs;
	if (_prev_vtxs) delete [] _prev_vtxs;
	if (_vtx_boxes) delete [] _vtx_boxes;

	if (_colors) delete [] _colors;
	if (_nrms) delete [] _nrms;

	if (_edges) delete [] _edges;
	if (_edg_boxes) delete [] _edg_boxes;

	if (_tris) delete [] _tris;
	if (_tri_edges) delete [] _tri_edges;
	if (_fac_boxes) delete [] _fac_boxes;
	
	if (_tri_nrms) delete [] _tri_nrms;
	if (_old_tri_nrms) delete [] _old_tri_nrms;
	if (_tri_flags) delete [] _tri_flags;

	if (_parts) delete [] _parts;
}

void
DeformModel::Display()
{
	glShadeModel(GL_SMOOTH);
	glEnableClientState( GL_VERTEX_ARRAY );
	glEnableClientState( GL_NORMAL_ARRAY );

	if (_colors) {
		glEnable(GL_COLOR_MATERIAL);
		glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);
		glEnableClientState( GL_COLOR_ARRAY );

		glColorPointer(3, GL_UNSIGNED_BYTE, 0, _colors);
	}

	glVertexPointer(3, GL_FLOAT, 0, _cur_vtxs);
	glNormalPointer(GL_FLOAT, 0, _nrms);

	glDrawElements( GL_TRIANGLES, _num_tri*3, GL_UNSIGNED_INT, _tris);

	glDisableClientState( GL_VERTEX_ARRAY );
	glDisableClientState( GL_NORMAL_ARRAY );
	glDisableClientState( GL_COLOR_ARRAY );
	
	glDisable(GL_COLOR_MATERIAL);
}

bool
DeformModel::Deform(float t, float circle)
{
	if (t >= circle)
		return false;

	t = t/circle*_num_frame;
	unsigned int prev_frame = ((unsigned int)t)%_num_frame;
	unsigned int next_frame = (prev_frame+1)%_num_frame;

	if (next_frame < prev_frame)
		return false;

	tm.startTiming(9);
	tm.startTiming(11);
	UpdateVert(prev_frame, next_frame, t-prev_frame);
	UpdateTriNorm();

	tm.endTiming(11);
	tm.endTiming(9);

	return true;
}

#define swap(a, b) {\
	vec3f *tmp = a;\
	a = b;\
	b = tmp;\
}

inline vec3f interp(const vec3f &p1, const vec3f &p2, float t)
{
	return p1*(1-t)+p2*t;
}


inline vec3f update(vec3f &v1, vec3f &v2, vec3f &v3)
{
	vec3f s = (v2-v1);
	return s.cross(v3-v1);
}

inline vec3f
update(tri3f &tri, vec3f *vtxs)
{
	vec3f &v1 = vtxs[tri.id0()];
	vec3f &v2 = vtxs[tri.id1()];
	vec3f &v3 = vtxs[tri.id2()];

	return update(v1, v2, v3);
}

void
DeformModel::UpdateVert(unsigned int prev, unsigned int next, float t)
{
	vec3f *prev_pt = _vtxs+prev*_num_vtx;
	vec3f *next_pt = _vtxs+next*_num_vtx;

	swap(_prev_vtxs, _cur_vtxs);

#pragma omp parallel for schedule(guided, 100)
	for (int i=0; i<_num_vtx; i++) {
		_cur_vtxs[i] = interp(prev_pt[i], next_pt[i], t);
	}
#pragma omp barrier
}

void
DeformModel::UpdateTriNorm()
{
#pragma omp parallel for schedule(guided, 100)
	for (int i=0; i<_num_vtx; i++) {
		_vtx_boxes[i] = BOX(_cur_vtxs[i]) + _prev_vtxs[i];
	}

#pragma omp parallel for schedule(guided, 200)
	for (int i=0; i<_num_edge; i++) {
		unsigned int id0 = _edges[i].vid(0);
		unsigned int id1 = _edges[i].vid(1);

		_edg_boxes[i] = BOX(_cur_vtxs[id0]) + _prev_vtxs[id0] + _cur_vtxs[id1] + _prev_vtxs[id1];
	}

#pragma omp parallel for schedule(guided, 100)
	for (int i=0; i<_num_tri; i++) {
		unsigned int id0 = _tris[i].id0();
		unsigned int id1 = _tris[i].id1();
		unsigned int id2 = _tris[i].id2();

		_fac_boxes[i] = BOX(_cur_vtxs[id0]) + _prev_vtxs[id0] +
								_cur_vtxs[id1] + _prev_vtxs[id1] +
								_cur_vtxs[id2] + _prev_vtxs[id2];
	}
}

void
DeformModel::BuildBVH(bool ccd)
{
	_tree = new UpperTree(this, ccd);
}

void
DeformModel::RebuildBVH(bool ccd)
{
	_tree->rebuild(ccd);
}

float DeformModel::RefitBVH(bool ccd)
{
	_tree->refit(ccd);
	return 0.f;
}

void DeformModel::DisplayBVH(int level, bool upper)
{
	glDisable(GL_LIGHTING);
	if (_tree)
		_tree->visulization(level, upper);
	glEnable(GL_LIGHTING);
}

void DeformModel::ResetCounter()
{
	// reset results
	for (int i=0; i<16; i++)
		_counts[i].init();
}

void DeformModel::UpdateCollide()
{
	if (_tree == NULL)
		return;

	non_adj_list.clear();

	tm.startTiming(12);
	front_list.propogate(this);
	//::sort(non_adj_list.begin(), non_adj_list.end());
	tm.endTiming(12);

	do_pairs();
}

void DeformModel::SelfCollide(bool ccd)
{
	if (_tree == NULL)
		return;

	non_adj_list.clear();

	front_list.clear();

	tm.startTiming(1);
	_tree->self_collide();
	tm.endTiming(1);

	do_pairs();
}

void
DeformModel::do_pairs()
{
	tm.startTiming(4);
	do_non_adj_pairs();
	tm.endTiming(4);

	tm.startTiming(2);
	do_orphans();
	tm.endTiming(2);
}

void
DeformModel::do_non_adj_pairs()
{
	unsigned int id1, id2;


	int length = non_adj_list.size();
#pragma omp parallel for
	for (int i=0; i<length; i++) {
		non_adj_list[i].get_param(id1, id2);
		test_feature_0(id1, id2);
	}
}

unsigned int
DeformModel::Covertex_F(unsigned int id1, unsigned int id2, unsigned int &st1, unsigned int &st2)
{
	unsigned int keeps[4];
	unsigned int num = 0;
	for (int i=0; i<3; i++)
		for (int j=0; j<3; j++) {
			if (_tris[id1].id(i) == _tris[id2].id(j)) {
				if (num < 2) {
					keeps[num*2] = i;
					keeps[num*2+1] = j;
				}
				num++;
			}
		}

	assert(num <= 3);
	if (num == 1) {
		st1 = keeps[0];
		st2 = keeps[1];
	} else
	if (num == 2) {
		for (int i=0; i<3; i++) {
			if (i != keeps[0] && i != keeps[2]) {
				st1 = i;
			}
			if (i != keeps[1] && i != keeps[3]) {
				st2 = i;
			}
		}
	} else {
		st1 = st2 = 0;
	}

	return num;
}

