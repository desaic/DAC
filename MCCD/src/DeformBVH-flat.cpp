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

#include <stdlib.h>
#include <assert.h>

#include "DeformBVH.h"
#include "DeformModel.h"
#include "aap.h"

#include "timing.h"
extern CBVHTimer tm;

static DeformBVHTree *s_tree;
static bool s_part = false;

#define ID(a) (s_part ? s_tree->_tri_idxes[(a)] : (a))

extern float middle_xyz(char xyz, const vec3f &p1, const vec3f &p2, const vec3f &p3);

DeformBVHTree::DeformBVHTree(DeformModel *mdl, bool ccd, unsigned int part)
{
	_part = part;
	s_part = (_part != -1);

	Construct(mdl, ccd);
}

static unsigned int *s_idx_buffer;
static DeformBVHNode *s_current_node;
static DeformModel *s_mdl;

void
DeformBVHTree::Construct(DeformModel *mdl, bool ccd)
{
	BOX total;
	int count;

	if (_part == -1) {
		for (unsigned int i=0; i<mdl->_num_vtx; i++) {
			total += mdl->_cur_vtxs[i];
			total += mdl->_prev_vtxs[i];
		}

		count = mdl->_num_tri;
	} else {
		count = 0;
		for (unsigned int i=0; i<mdl->_num_tri; i++) {
			if (mdl->_parts[i] != _part)
				continue;
			else
				count++;

			for (int j=0; j<3; j++) {
				total += mdl->_cur_vtxs[mdl->_tris[i].id(j)];
				total += mdl->_prev_vtxs[mdl->_tris[i].id(j)];
			}
		}
	}

	assert(mdl->_tri_boxes == NULL);
	assert(mdl->_tri_centers == NULL);

	mdl->_tri_boxes = new BOX[count];
	mdl->_tri_centers = new vec3f[count];

	if (_part != -1)
		_tri_idxes = new unsigned int[count];
	else
		_tri_idxes = NULL;

	aap pln(total);

	s_idx_buffer = new unsigned int[count];
	unsigned int left_idx = 0, right_idx = count, tri_idx = 0;

	for (unsigned int i=0; i<mdl->_num_tri; i++) {
		if (_part != -1 && mdl->_parts[i] != _part)
			continue;

		tri_idx++;
		if (_part != -1)
			_tri_idxes[tri_idx-1] = i;

		vec3f &p1 = mdl->_cur_vtxs[mdl->_tris[i].id0()];
		vec3f &p2 = mdl->_cur_vtxs[mdl->_tris[i].id1()];
		vec3f &p3 = mdl->_cur_vtxs[mdl->_tris[i].id2()];
		vec3f &pp1 = mdl->_prev_vtxs[mdl->_tris[i].id0()];
		vec3f &pp2 = mdl->_prev_vtxs[mdl->_tris[i].id1()];
		vec3f &pp3 = mdl->_prev_vtxs[mdl->_tris[i].id2()];
		
		mdl->_tri_centers[tri_idx-1].set_value(
				(middle_xyz(0, p1, p2, p3)+middle_xyz(0, pp1, pp2, pp3))*0.5f,
				(middle_xyz(1, p1, p2, p3)+middle_xyz(1, pp1, pp2, pp3))*0.5f,
				(middle_xyz(2, p1, p2, p3)+middle_xyz(2, pp1, pp2, pp3))*0.5f);

		if (pln.inside(mdl->_tri_centers[tri_idx-1]))
			s_idx_buffer[left_idx++] = tri_idx-1;
		else
			s_idx_buffer[--right_idx] = tri_idx-1;

		mdl->_tri_boxes[tri_idx-1] += p1;
		mdl->_tri_boxes[tri_idx-1] += p2;
		mdl->_tri_boxes[tri_idx-1] += p3;
		mdl->_tri_boxes[tri_idx-1] += pp1;
		mdl->_tri_boxes[tri_idx-1] += pp2;
		mdl->_tri_boxes[tri_idx-1] += pp3;
	
	}

	_nodes = new DeformBVHNode[count*2-1];
	_nodes[0]._box = total;
	_nodes[0]._count = count;
	s_current_node = _nodes+3;
	s_tree = this;

	if (count == 1) {
		_nodes[0]._child = ID(0);
	} else {
		_nodes[0]._child = -1;

		if (left_idx == 0 || left_idx == count)
			left_idx = count/2;

		s_mdl = mdl;
		_nodes[0].getLeftChild()->Construct(s_idx_buffer, left_idx);
		_nodes[0].getRightChild()->Construct(s_idx_buffer+left_idx, count-left_idx);
	}

	_mdl = mdl;

	delete [] mdl->_tri_boxes;
	mdl->_tri_boxes = NULL;
	delete [] mdl->_tri_centers;
	mdl->_tri_centers = NULL;
	delete [] s_idx_buffer;
	s_idx_buffer = NULL;
}

DeformBVHTree::~DeformBVHTree()
{
	delete [] _nodes;
}

//#################################################################
// called by root
DeformBVHNode::DeformBVHNode()
{
	_child = 	_count = 0;
}

DeformBVHNode::~DeformBVHNode()
{
	NULL;
}

void
DeformBVHNode::Construct(unsigned int id)
{
	_count = 1;
	_child = id;
	//_box = s_mdl->_tri_boxes[id];
}

void
DeformBVHNode::Construct(unsigned int *lst, unsigned int lst_num)
{
	_count = lst_num;
	for (unsigned int t=0; t<lst_num; t++)
		_box += s_mdl->_tri_boxes[lst[t]];

	if (lst_num == 1) {
		_child = ID(lst[0]);
		_box = s_mdl->_tri_boxes[lst[0]];
	} else { // try to split
		_child = this-s_current_node;
		s_current_node += 2;

		if (lst_num == 2) {// must split
			getLeftChild()->Construct(ID(lst[0]));
			getRightChild()->Construct(ID(lst[1]));
		} else {
			aap pln(_box);
			unsigned int left_idx = 0, right_idx = lst_num-1;
			for (unsigned int t=0; t<lst_num; t++) {
				int i=lst[left_idx];
				if (pln.inside(s_mdl->_tri_centers[i]))
					left_idx++;
				else {// swap it
					unsigned int tmp = i;
					lst[left_idx] = lst[right_idx];
					lst[right_idx--] = tmp;
				}
			}

			int half = lst_num/2;

			if (left_idx == 0 || left_idx == lst_num) {
				getLeftChild()->Construct(lst, half);
				getRightChild()->Construct(lst+half, lst_num-half);
			} else {
				getLeftChild()->Construct(lst, left_idx);
				getRightChild()->Construct(lst+left_idx, lst_num-left_idx);
			}
		}
	}
}
