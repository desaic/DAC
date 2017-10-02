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

#include "DeformModel.h"
#include "plyfile.h"
#pragma warning(disable: 4996)

#include <vector>
#include <list>
#include <algorithm>

typedef struct PLYVertex{
	float coords[3];
	unsigned char color[3];
	void *other_props;		
} PLYVertex;

typedef struct PLYFace{
	unsigned char nverts;
	int *verts;		
	void *other_props;
} PLYFace;

#include "logger.h"
extern Logger *logger;

PlyProperty vert_props[] = { /* list of property information for a vertex */
	{"x", PLY_FLOAT, PLY_FLOAT, 0, 0, 0, 0, 0},
	{"y", PLY_FLOAT, PLY_FLOAT, 4, 0, 0, 0, 0},
	{"z", PLY_FLOAT, PLY_FLOAT, 8, 0, 0, 0, 0},
	{"red", PLY_UCHAR, PLY_UCHAR, (int)offsetof(PLYVertex,color[0]), 0, 0, 0, 0},
	{"green", PLY_UCHAR, PLY_UCHAR, (int)offsetof(PLYVertex,color[1]), 0, 0, 0, 0},
	{"blue", PLY_UCHAR, PLY_UCHAR, (int)offsetof(PLYVertex,color[2]), 0, 0, 0, 0},
};

PlyProperty face_props[] = { /* list of property information for a vertex */
	{"vertex_indices", PLY_INT, PLY_INT, offsetof(PLYFace,verts), 1, PLY_UCHAR, PLY_UCHAR, offsetof(PLYFace,nverts)},
};

#define swapI(a, b) {\
	unsigned int tmp = a;\
	a = b;\
	b = tmp;\
}

char
DeformModel::get_status_1(unsigned int id1, unsigned int id2, unsigned int st1, unsigned int st2)
{
	unsigned int e0 = _tri_edges[id1].id(st1);
	unsigned int e00 = _tri_edges[id1].id((st1+2)%3);

	unsigned int fid1 = _edges[e0].fid(0);
	unsigned int f1_mates = (fid1 == id1) ? _edges[e0].fid(1) : fid1;

	unsigned int e1 = _tri_edges[id2].id(st2);
	unsigned int e11 = _tri_edges[id2].id((st2+2)%3);
	if (Covertex_E(e0, e1)) {
		swapI(e1, e11);
	}
	unsigned int fid2 = _edges[e11].fid(0);
	unsigned int f2_mates = (fid2 == id2) ? _edges[e11].fid(1) : fid2;

	if (f1_mates != -1)
		if (f2_mates != -1)
			return 3; // just skip it
		else
			return 2; //group 2 only
	else
		if (f2_mates != -1)
			return 1; //group 1 only
		else
			return 0; //both group
}

char
DeformModel::get_status_2(unsigned int id1, unsigned int id2, unsigned int st1, unsigned int st2)
{
	unsigned int edg1 = _tri_edges[id1].id((st1+1)%3);
	unsigned int fid1 = _edges[edg1].fid(0);
	unsigned int f1_mates = (fid1 == id1) ? _edges[edg1].fid(1) : fid1;
		
	unsigned int edg2 = _tri_edges[id2].id((st2+1)%3);
	unsigned int fid2 = _edges[edg2].fid(0);
	unsigned int f2_mates = (fid2 == id2) ? _edges[edg2].fid(1) : fid2;


	if (f1_mates != -1)
		if (f2_mates != -1)
			return 3; //just skip it ...
		else
			return 2; //group 2 only
	else
		if (f2_mates != -1)
			return 1; //group 1 only
		else
			return 0; //both group
}

void
DeformModel::BufferAdjacent()
{
	adj_pair_list adj_2_list, adj_1_list;

	adj_2_list.clear();
	for (unsigned int i=0; i<_num_edge; i++) {
		unsigned int id1 = _edges[i].fid(0);
		unsigned int id2 = _edges[i].fid(1);

		if (id1 == -1 || id2 == -1) continue;
		if (id1 < id2) swapI(id1, id2);

		unsigned int st1 = 0, st2 = 0;
		unsigned int cov = Covertex_F(id1, id2, st1, st2);
		//20090903 ...
		//assert(cov == 2);

		char status = get_status_2(id1, id2, st1, st2);
		//if (status == 3) continue;

		adj_2_list.push_back(adjacent_pair(id1, id2, st1, st2, status));
	}

	adj_1_list.clear();
	for (unsigned int i=0; i<_num_vtx; i++) {
		for (id_list::iterator it1=_vtx_fids[i].begin(); it1!=_vtx_fids[i].end(); it1++) {
			for (id_list::iterator it2=it1; it2!=_vtx_fids[i].end(); it2++) {
				if (it1 == it2) continue;

				unsigned int id1 = *it1;
				unsigned int id2 = *it2;
				if (id1 < id2) swapI(id1, id2);

				unsigned int st1 = 0, st2 = 0;
				unsigned int cov = Covertex_F(id1, id2, st1, st2);
				if (cov == 2) continue;

				char status = get_status_1(id1, id2, st1, st2);
				//if (status == 3) continue;

				adj_1_list.push_back(adjacent_pair(id1, id2, st1, st2, status));
			}
		}
	}

	sort(adj_1_list.begin(), adj_1_list.end());
	adj_1_list.erase(unique(adj_1_list.begin(), adj_1_list.end()), adj_1_list.end());

	get_orphans(adj_2_list, adj_1_list);
}

DeformModel::DeformModel(double * verts, int nv, int * tidx, int nt,
  unsigned int num_frame)
{
  Init();
  
  _num_vtx = nv;
  _vtxs = new vec3f[_num_vtx*_num_frame];
  _cur_vtxs = new vec3f[_num_vtx];
  _prev_vtxs = new vec3f[_num_vtx];
  _vtx_boxes = new BOX[_num_vtx];
  _nrms = new vec3f[_num_vtx];
  _vtx_fids = new id_list[_num_vtx];
  std::cout << "Vtx # = " << _num_vtx << endl;
  
  _num_frame = num_frame;
  for (int j = 0; j < nv; j++) {
    _vtxs[j] = vec3f(verts[3 * j], verts[3 * j + 1], verts[3 * j + 2]);
    _prev_vtxs[j] = _cur_vtxs[j] = _vtxs[j];
  }


  list<edge2f> edgelist_temp;
  vector<tri3f> trilist_temp;

  for (int j = 0; j < nt; j++) {

    int id0 = tidx[3 * j];
    int id1 = tidx[3 * j + 1];
    int id2 = tidx[3 * j + 2];
    tri3f tri(id0, id1, id2);

    // insert triangle into list
    trilist_temp.push_back(tri);
    unsigned int fid = (unsigned int)trilist_temp.size() - 1;

    edgelist_temp.push_back(edge2f(id0, id1, fid));
    edgelist_temp.push_back(edge2f(id1, id2, fid));
    edgelist_temp.push_back(edge2f(id2, id0, fid));
  }
  edgelist_temp.sort();

  list<edge2f> edge_unqie;
  for (list<edge2f>::iterator it = edgelist_temp.begin(); it != edgelist_temp.end(); it++) {
    if (!edge_unqie.empty() && *it == edge_unqie.back()) { // find duplicated with other fid
      unsigned int fid = (*it).fid(0);
      assert(fid != -1);
      edge_unqie.back().set_fid2(fid);
    }
    else
      edge_unqie.push_back(*it);
  }

  edgelist_temp.clear();
  vector<edge2f> edge_array;

  _num_edge = (unsigned int)edge_unqie.size();
  _edges = new edge2f[_num_edge];
  _edg_boxes = new BOX[_num_edge];

  unsigned int t = 0;
  for (list<edge2f>::iterator it = edge_unqie.begin(); it != edge_unqie.end(); it++)
  {
    _edges[t++] = *it;
    edge_array.push_back(*it);
  }

  // copy over temp list to static array
  _num_tri = (unsigned int)trilist_temp.size();
  std::cout << "Allocating " << _num_tri * sizeof(tri3f) << " bytes of storage for triangles." << endl;
  _tris = new tri3f[_num_tri];
  _tri_nrms = new vec3f[_num_tri];
  _old_tri_nrms = NULL;

  _tri_edges = new tri3e[_num_tri];
  _fac_boxes = new BOX[_num_tri];
  _tri_flags = new char[_num_tri];

  vector <edge2f>::iterator first = edge_array.begin();
  vector <edge2f>::iterator last = edge_array.end();
  for (t = 0; t < _num_tri; t++) {
    _tris[t] = trilist_temp[t];

    vector <edge2f>::iterator it1 = lower_bound(first, last, edge2f(_tris[t].id0(), _tris[t].id1(), 0));
    vector <edge2f>::iterator it2 = lower_bound(first, last, edge2f(_tris[t].id1(), _tris[t].id2(), 0));
    vector <edge2f>::iterator it3 = lower_bound(first, last, edge2f(_tris[t].id2(), _tris[t].id0(), 0));

    _tri_edges[t].set(it1 - first, it2 - first, it3 - first);
  }

  UpdateTriNorm();

  for (unsigned t = 0; t < _num_tri; t++)
    for (int i = 0; i < 3; i++) {
      unsigned int vid = _tris[t].id(i);
      if (vid > _num_vtx) {
        std::cout << "debug\n";
      }
      _vtx_fids[vid].push_back(t);
    }

  BufferAdjacent();
}

DeformModel::DeformModel(char *fname, unsigned int num_frame, float ply_scale)
{
	Init();

	char ply_fname[256];
	
	_num_frame = num_frame;
	for (unsigned int cur_f = 0; cur_f < _num_frame; cur_f++) {
		sprintf(ply_fname, "%s%d.ply", fname, cur_f);
		FILE *fp = fopen(ply_fname, "rb");
		assert(fp);

		// PLY object:
		PlyFile *ply;

		// PLY properties:
		char **elist;
		int nelems;

		// hand over the stream to the ply functions:
		ply = ply_read(fp, &nelems, &elist);
		assert(ply);

		int file_type;
		float version;		
		ply_get_info(ply, &version, &file_type);

		for (int i=0; i<nelems; i++) {
			char *elem_name = elist[i];

			int num_elems, nprops;
			PlyProperty **plist = ply_get_element_description(ply, elem_name, &num_elems, &nprops);

			bool has_vertex_x = false, has_vertex_y = false, has_vertex_z = false, has_colors = false;
			unsigned char color_components = 0;

			// this is a vertex:
			if (equal_strings ("vertex", elem_name)) {
				for (int j=0; j<nprops; j++)
				{
					if (equal_strings("x", plist[j]->name))
					{
						ply_get_property (ply, elem_name, &vert_props[0]);  /* x */
						has_vertex_x = TRUE;
					}
					else if (equal_strings("y", plist[j]->name))
					{
						ply_get_property (ply, elem_name, &vert_props[1]);  /* y */
						has_vertex_y = TRUE;
					}
					else if (equal_strings("z", plist[j]->name))
					{
						ply_get_property (ply, elem_name, &vert_props[2]);  /* z */
						has_vertex_z = TRUE;
					}
					else if (equal_strings("red", plist[j]->name))
					{
						ply_get_property (ply, elem_name, &vert_props[3]);  /* z */
						color_components++;
					}
					else if (equal_strings("green", plist[j]->name))
					{
						ply_get_property (ply, elem_name, &vert_props[4]);  /* z */
						color_components++;
					}
					else if (equal_strings("blue", plist[j]->name))
					{
						ply_get_property (ply, elem_name, &vert_props[5]);  /* z */
						color_components++;
					}
				}

				has_colors = color_components == 3;
				// test for necessary properties
				if ((!has_vertex_x) || (!has_vertex_y) || (!has_vertex_z))
				{
					logger->update("Warning: Vertex with less than 3 coordinated detected. Output will most likely be corrupt!");
					continue;
				}

				// must be first frame, initialize structures:
				if (_num_vtx == 0) {
					_num_vtx = num_elems;
					_vtxs = new vec3f[_num_vtx*_num_frame];
					_cur_vtxs = new vec3f[_num_vtx];
					_prev_vtxs = new vec3f[_num_vtx];
					_vtx_boxes = new BOX[_num_vtx];
					_nrms = new vec3f[_num_vtx];
					_vtx_fids = new id_list[_num_vtx];
					cout << "Vtx # = " << _num_vtx << endl;

					if (has_colors)
						_colors = new color3[_num_vtx];
				}

				// grab all the vertex elements
				PLYVertex plyNewVertex;		
				for (int j=0; j<num_elems; j++) {
					ply_get_element(ply, (void *)&plyNewVertex);								

					if (has_colors && cur_f == 0) {
						_colors[j].set(plyNewVertex.color);
					}

					_vtxs[cur_f*_num_vtx+j] = vec3f(plyNewVertex.coords) * ply_scale;
					if (cur_f == 0) {
						_prev_vtxs[j] = _cur_vtxs[j] = _vtxs[j];
					}
					
					if (j != 0 && j%1000000 == 0) {
						cout << " - " << j << " of " << num_elems << " loaded." << endl;					
					}				
				}
			}
			
			// this is a face (and, hopefully, a triangle):
			else if (equal_strings ("face", elem_name) && _tris == NULL) {
				// I need this for..., otherwise error ...
				for (int j=0; j<nprops; j++)
				{
					if (equal_strings("vertex_indices", plist[j]->name))
					{
						ply_get_property (ply, elem_name, &face_props[0]);  /* vertex_indices */
					}
				}

				/* grab all the face elements */
				PLYFace plyFace;	
				plyFace.other_props = NULL;			
								
				list<edge2f> edgelist_temp;
				vector<tri3f> trilist_temp;

				for (int j = 0; j < num_elems; j++) {
					ply_get_element(ply, (void *)&plyFace);
					for (int fi = 0; fi < plyFace.nverts-2; fi++) {
						//
						// make a triangle in our format from PLY face + vertices
						//						
						// copy vertex indices
						unsigned int id0, id1, id2;

						id0 = plyFace.verts[0];
						id1 = plyFace.verts[fi+1];
						id2 = plyFace.verts[fi+2];

						tri3f tri(id0, id1, id2);

						// insert triangle into list
						trilist_temp.push_back(tri);
						unsigned int fid = (unsigned int)trilist_temp.size()-1;

						edgelist_temp.push_back(edge2f(id0, id1, fid));
						edgelist_temp.push_back(edge2f(id1, id2, fid));
						edgelist_temp.push_back(edge2f(id2, id0, fid));
					}
					free(plyFace.verts);												

					if (j != 0 && j%500000 == 0) {
						cout << " - " << j << " of " << num_elems << " loaded." << endl;					
					}
				}

				edgelist_temp.sort();
				
				list<edge2f> edge_unqie;
				for (list<edge2f>::iterator it=edgelist_temp.begin(); it!=edgelist_temp.end(); it++) {
					if (!edge_unqie.empty() && *it == edge_unqie.back()) { // find duplicated with other fid
							unsigned int fid = (*it).fid(0);
							assert(fid != -1);
							edge_unqie.back().set_fid2(fid);
						} else
							edge_unqie.push_back(*it);
				}

				edgelist_temp.clear();
				vector<edge2f> edge_array;

				_num_edge = (unsigned int)edge_unqie.size();
				_edges = new edge2f[_num_edge];
				_edg_boxes = new BOX[_num_edge];

				unsigned int t=0;
				for (list<edge2f>::iterator it=edge_unqie.begin(); it != edge_unqie.end(); it++)
				{
					_edges[t++] = *it;
					edge_array.push_back(*it);
				}

				// copy over temp list to static array
				_num_tri = (unsigned int)trilist_temp.size();
				cout << "Allocating " << _num_tri*sizeof(tri3f) << " bytes of storage for triangles." << endl;			
				_tris = new tri3f[_num_tri];
				_tri_nrms = new vec3f[_num_tri];
				_old_tri_nrms = NULL;

				_tri_edges = new tri3e[_num_tri];
				_fac_boxes = new BOX[_num_tri];
				_tri_flags = new char[_num_tri];

				vector <edge2f>::iterator first = edge_array.begin();
				vector <edge2f>::iterator last = edge_array.end();
				for (t = 0; t < _num_tri; t++) {
					_tris[t] = trilist_temp[t];

					vector <edge2f>::iterator it1 = lower_bound(first, last, edge2f(_tris[t].id0(), _tris[t].id1(), 0));
					vector <edge2f>::iterator it2 = lower_bound(first, last, edge2f(_tris[t].id1(), _tris[t].id2(), 0));
					vector <edge2f>::iterator it3 = lower_bound(first, last, edge2f(_tris[t].id2(), _tris[t].id0(), 0));

					_tri_edges[t].set(it1-first, it2-first, it3-first);
				}

				cout << "Edge # = " << _num_edge << endl;
				cout << "Tri # = " << _num_tri << endl;
			}

			else // otherwise: skip all further
				NULL;
		}

		// PLY parsing ended, clean up vertex buffer and close the file		
		ply_close(ply);
		fclose(fp);
	}

	UpdateTriNorm();

	for (unsigned t = 0; t < _num_tri; t++)
		for (int i=0; i<3; i++) {
			unsigned int vid = _tris[t].id(i);
			_vtx_fids[vid].push_back(t);
		}

	BufferAdjacent();
}
