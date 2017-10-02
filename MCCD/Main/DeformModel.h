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

#pragma once
#include <stdlib.h>
#include <omp.h>

#include "vec3f.h"

#include <vector>
#include <set>

using namespace std;
typedef vector<unsigned int> id_list;

#include "tri_pair.h"
typedef vector<adjacent_pair> adj_pair_list;

class edge2f {
	unsigned int _vids[2];
	unsigned int _fids[2];

	FORCEINLINE void set(unsigned int id0, unsigned int id1) {
		if (id0 > id1) {
			_vids[0] = id0;
			_vids[1] = id1;
		} else {
			_vids[1] = id0;
			_vids[0] = id1;
		}
	}
public:
	FORCEINLINE edge2f() {
		_vids[0] = _vids[1] = UINT_MAX;
		_fids[0] = _fids[1] = UINT_MAX;
	}

	FORCEINLINE edge2f(unsigned id0, unsigned int id1, unsigned int fid) {
		set(id0, id1);
		_fids[0] = fid;
		_fids[1] = UINT_MAX;
	}

	FORCEINLINE void set_fid2(unsigned id) {
		_fids[1] = id;
	}

	FORCEINLINE unsigned int vid(int i) { return _vids[i]; }
	FORCEINLINE unsigned int fid(int i) { return _fids[i]; }

	FORCEINLINE bool operator == (const edge2f &other) const {
		return (_vids[0] == other._vids[0] && _vids[1] == other._vids[1]);
	}

	FORCEINLINE bool operator < (const edge2f &other) const {
		if (_vids[0] == other._vids[0])
			return _vids[1] < other._vids[1];
		else
			return _vids[0] < other._vids[0];
	}
};

class tri3e {
	unsigned int _ids[3];

public:
	FORCEINLINE tri3e() {
		_ids[0] = _ids[1] = _ids[2] = UINT_MAX;
	}

	FORCEINLINE tri3e(unsigned int id0, unsigned int id1, unsigned id2) {
		set(id0, id1, id2);
	}

	FORCEINLINE void set(unsigned int id0, unsigned int id1, unsigned int id2) {
		_ids[0] = id0;
		_ids[1] = id1;
		_ids[2] = id2;
	}

	FORCEINLINE unsigned int id(int i) { return _ids[i]; }
};

class tri3f {
	unsigned int _ids[3];

public:
	FORCEINLINE tri3f() {
		_ids[0] = _ids[1] = _ids[2] = UINT_MAX;
	}

	FORCEINLINE tri3f(unsigned int id0, unsigned int id1, unsigned id2) {
		set(id0, id1, id2);
	}

	FORCEINLINE void set(unsigned int id0, unsigned int id1, unsigned int id2) {
		_ids[0] = id0;
		_ids[1] = id1;
		_ids[2] = id2;
	}

	FORCEINLINE unsigned int id(int i) { return _ids[i]; }
	FORCEINLINE unsigned int id0() {return _ids[0];}
	FORCEINLINE unsigned int id1() {return _ids[1];}
	FORCEINLINE unsigned int id2() {return _ids[2];}
};

class color3 {
public:
	unsigned char _rgbs[3];

	FORCEINLINE color3() {
		_rgbs[0] = _rgbs[1] = _rgbs[2] = 0;
	}

	FORCEINLINE void set(unsigned char rgb[]) {
		_rgbs[0] = rgb[0];
		_rgbs[1] = rgb[1];
		_rgbs[2] = rgb[2];
	}

	FORCEINLINE void set(unsigned char r, unsigned char g, unsigned char b) {
		_rgbs[0] = r;
		_rgbs[1] = g;
		_rgbs[2] = b;
	}

};

class UpperTree;

#include "box.h"

class DeformModel {
	unsigned int _num_vtx;
	unsigned int _num_frame;

	vec3f *_vtxs;
	vec3f *_cur_vtxs;
	vec3f *_prev_vtxs;
	id_list *_vtx_fids;
	BOX *_vtx_boxes;

	color3 *_colors;
	vec3f *_nrms;

	unsigned int _num_tri;
	tri3f *_tris;

	tri3e *_tri_edges;
	BOX *_fac_boxes;

	vec3f *_tri_nrms;
	vec3f *_old_tri_nrms;

	char *_tri_flags;

	unsigned int _num_edge;
	edge2f *_edges;
	BOX *_edg_boxes;

	// for building BVH
	UpperTree *_tree;

	vec3f *_tri_centers;
	BOX *_tri_boxes;

	// for collide
	struct {
	unsigned int _num_box_tests;
	unsigned int _num_vf_test;
	unsigned int _num_ee_test;
	unsigned int _num_vf_true;
	unsigned int _num_ee_true;

	void init() {
		_num_box_tests =
		_num_vf_test =
		_num_ee_test =
		_num_vf_true =
		_num_ee_true = 0;

	}
	} _counts[16];

	unsigned int _num_parts;
	unsigned int *_parts;

	void do_pairs();

	void do_non_adj_pairs();

	void BufferAdjacent();
	char get_status_1(unsigned int id1, unsigned int id2, unsigned int st1, unsigned int st2);
	char get_status_2(unsigned int id1, unsigned int id2, unsigned int st1, unsigned int st2);

public:
	DeformModel();
	DeformModel(char *fname, unsigned int num_frame, float ply_scale);
	DeformModel(char *fname, unsigned int num_frame, unsigned int num_bodys, float ply_scale);
	DeformModel(char *fname, unsigned int num_frame); // aff
  DeformModel(double * verts, int nv, int * tidx, int nt,
    unsigned int num_frame);
  ~DeformModel();

	unsigned int GetFrames() { return _num_frame; }
	void Display();
	bool Deform(float t, float circle);

	void Init();
	void Clean();

	void UpdateVert(unsigned int prev, unsigned int next, float t);
	void UpdateVtxNorm();
	void UpdateTriNorm();


	void BuildBVH(bool ccd);
	void RebuildBVH(bool ccd);
	float RefitBVH(bool ccd);
	void DisplayBVH(int, bool upper);

	void ResetCounter();
	void SelfCollide(bool ccd);

	void UpdateCollide();

	void Decompose();

	FORCEINLINE int NumTri() { return _num_tri; }

	FORCEINLINE int NumBoxTest()
	{
		int sum = 0;
		for (int i=0; i<16; i++)
		sum += _counts[i]._num_box_tests;
		return sum;
	}

	FORCEINLINE int NumVFTest()
	{
		int sum = 0;
		for (int i=0; i<16; i++)
		sum += _counts[i]._num_vf_test;
		return sum;
	}

	FORCEINLINE int NumEETest()
	{
		int sum = 0;
		for (int i=0; i<16; i++)
		sum += _counts[i]._num_ee_test;
		return sum;
	}

	FORCEINLINE int NumVFTrue()
	{
		int sum = 0;
		for (int i=0; i<16; i++)
		sum += _counts[i]._num_vf_true;
		return sum;
	}

	FORCEINLINE int NumEETrue()
	{
		int sum = 0;
		for (int i=0; i<16; i++)
		sum += _counts[i]._num_ee_true;
		return sum;
	}


	FORCEINLINE bool Covertex_E(unsigned int e1, unsigned int e2) {
		for (int i=0; i<2; i++)
			for (int j=0; j<2; j++)
				if (_edges[e1].vid(i) == _edges[e2].vid(j))
					return true;

		return false;
	}

	FORCEINLINE bool Coedge_F(unsigned int id1, unsigned int id2) {
		unsigned int num = 0;
		for (int i=0; i<3; i++)
			for (int j=0; j<3; j++) {
				if (_tris[id1].id(i) == _tris[id2].id(j)) {
					num++;
				}
			}

		return num > 1;
	}

	FORCEINLINE bool Covertex_F(unsigned int id1, unsigned int id2) {
		for (int i=0; i<3; i++)
			for (int j=0; j<3; j++) {
				if (_tris[id1].id(i) == _tris[id2].id(j))
					return true;
			}

		return false;
	}

	unsigned int Covertex_F(unsigned int id1, unsigned int id2, unsigned int &st1, unsigned int &st2);

	friend class DeformBVHTree;
	friend class DeformBVHNode;
	friend class UpperTree;
	friend class UpperNode;
	friend class ee_vf_manager;
	friend class CollideBody;
	friend class bvh_front_node;

	float intersect_vf(unsigned int fid1, unsigned int vid2, unsigned int fid2);
	float intersect_ee(unsigned int e1, unsigned int e2, unsigned int fid1, unsigned int fid2);

	float intersect_vf(unsigned int fid, unsigned int vid);
	float intersect_ee(unsigned int e1, unsigned int e2);
	float do_vf(unsigned int fid, unsigned int vid);
	float do_ee(unsigned int e1, unsigned int e2);

	void test_feature_0(unsigned id1, unsigned int id2);

	void do_orphans();
	void load_orphans();
	void get_orphans(adj_pair_list& adj_2_list, adj_pair_list& adj_1_list);
	void get_feature_1(unsigned int, unsigned int, unsigned int, unsigned int);
	void get_feature_2(unsigned int, unsigned int, unsigned int, unsigned int);
	void insert_ee(unsigned int, unsigned int);
	void insert_vf(unsigned int, unsigned int);
	bool test_orphan_ee(unsigned int e1, unsigned int e2);
	bool test_orphan_vf(unsigned int f, unsigned int v);
};