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
#include "timing.h"
extern CBVHTimer tm;

static DeformModel *s_mdl1, *s_mdl2;

#include "tri_pair.h"
extern non_adjacent_pair_list non_adj_list;

#include "bvh_front.h"
extern bvh_front_list front_list;

void
DeformBVHNode::append_collide(DeformBVHNode *other, bvh_front_list &appended, non_adjacent_pair_list &result, int tid)
{
	if (isLeaf() && other->isLeaf()) {
		bool cov = s_mdl1->Covertex_F(getTriID(), other->getTriID());

		if (!cov) {
			appended.push_back(bvh_front_node(this, other));

			s_mdl1->_counts[tid]._num_box_tests++;

			if (!_box.overlaps(other->_box))
				return;

			result.push_back(non_adjacent_pair(getTriID(), other->getTriID()));
		}

		return;
	}

	s_mdl1->_counts[tid]._num_box_tests++;

	if (!_box.overlaps(other->_box)) {
		appended.push_back(bvh_front_node(this, other));
		return;
	}

	if (isLeaf()) {
		append_collide(other->getLeftChild(), appended, result, tid);
		append_collide(other->getRightChild(), appended, result, tid);
	} else {
		getLeftChild()->append_collide(other, appended, result, tid);
		getRightChild()->append_collide(other, appended, result, tid);
	}
}


void
DeformBVHTree::visulization(int level, bool dummy)
{
	getRoot()->visulization(level);
}

static DeformModel *s_mdl;

void
DeformBVHNode::getChildren(DeformBVHNode *&n1, DeformBVHNode *&n2, DeformBVHNode *&n3, DeformBVHNode *&n4)
{
	n1 = getLeftChild()->getLeftChild();
	n2 = getLeftChild()->getRightChild();
	n3 = getRightChild()->getLeftChild();
	n4 = getRightChild()->getRightChild();
}

void
DeformBVHNode::mergeBox(DeformBVHNode *n1, DeformBVHNode *n2, DeformBVHNode *n3, DeformBVHNode *n4)
{
	getLeftChild()->_box = n1->_box + n2->_box;
	getRightChild()->_box = n3->_box + n4->_box;
	_box = getLeftChild()->_box + getRightChild()->_box;
}

float
DeformBVHTree::refit(bool openmp)
{
	s_mdl = _mdl;

	if (openmp == false) {
		getRoot()->refit();
		return 0.f;
	}

	DeformBVHNode *n1, *n2, *n3, *n4;
	getRoot()->getChildren(n1, n2, n3, n4);
	DeformBVHNode *n11, *n12, *n13, *n14;
	n1->getChildren(n11, n12, n13, n14);
	DeformBVHNode *n21, *n22, *n23, *n24;
	n2->getChildren(n21, n22, n23, n24);
	DeformBVHNode *n31, *n32, *n33, *n34;
	n3->getChildren(n31, n32, n33, n34);
	DeformBVHNode *n41, *n42, *n43, *n44;
	n4->getChildren(n41, n42, n43, n44);

#pragma omp parallel sections
		{
 #pragma omp section
		n11->refit();
 #pragma omp section
		n12->refit();
 #pragma omp section
		n13->refit();
 #pragma omp section
		n14->refit();
 #pragma omp section
		n21->refit();
 #pragma omp section
		n22->refit();
 #pragma omp section
		n23->refit();
 #pragma omp section
		n24->refit();
 #pragma omp section
		n31->refit();
 #pragma omp section
		n32->refit();
 #pragma omp section
		n33->refit();
 #pragma omp section
		n34->refit();
 #pragma omp section
		n41->refit();
 #pragma omp section
		n42->refit();
 #pragma omp section
		n43->refit();
 #pragma omp section
		n44->refit();
		}

	n1->mergeBox(n11, n12, n13, n14);
	n2->mergeBox(n21, n22, n23, n24);
	n3->mergeBox(n31, n32, n33, n34);
	n4->mergeBox(n41, n42, n43, n44);
	getRoot()->mergeBox(n1, n2, n3, n4);

	return 0.f;
}

void
DeformBVHTree::collide(DeformBVHTree *other)
{
	s_mdl1 = _mdl;
	s_mdl2 = other->_mdl;

	getRoot()->collide(other->getRoot());
}

void
DeformBVHTree::self_collide()
{
	s_mdl1 = _mdl;
	s_mdl2 = _mdl;

	getRoot()->self_collide();
	return;
}

BOX
DeformBVHTree::box()
{
	return getRoot()->_box;
}

void
DeformBVHNode::visulization(int level)
{
	if (isLeaf())
		_box.visulization();
	else
		if ((level > 0)) {
			if (level == 1)
				_box.visulization();
			else {
				if (getLeftChild())
					getLeftChild()->visulization(level-1);
				if (getRightChild())
					getRightChild()->visulization(level-1);
			}
		}
}

inline vec3f norm(vec3f &p1, vec3f &p2, vec3f &p3)
{
	vec3f s = p2-p1;
	vec3f t = p3-p1;
	vec3f n = s.cross(t);
	return n;
}

void
DeformBVHNode::refit()
{
	if (isLeaf()) {
		_box = s_mdl->_fac_boxes[getTriID()];
	} else {
		getLeftChild()->refit();
		getRightChild()->refit();

		_box = getLeftChild()->_box + getRightChild()->_box;
	}
}

void
DeformBVHNode::self_collide()
{
	if (isLeaf())
		return;

	getLeftChild()->self_collide();
	getRightChild()->self_collide();
	getLeftChild()->collide(getRightChild());
}


void
DeformBVHNode::collide(DeformBVHNode *other)
{
	if (isLeaf() && other->isLeaf()) {
		bool cov = s_mdl1->Covertex_F(getTriID(), other->getTriID());

		if (!cov) {
			front_list.push_back(bvh_front_node(this, other));

			s_mdl1->_counts[0]._num_box_tests++;
			if (!_box.overlaps(other->_box))
				return;

			non_adj_list.push_back(non_adjacent_pair(getTriID(), other->getTriID()));
		}

		return;
	}

	s_mdl1->_counts[0]._num_box_tests++;
	if (!_box.overlaps(other->_box)) {
		front_list.push_back(bvh_front_node(this, other));
		return;
	}

	if (isLeaf()) {
		collide(other->getLeftChild());
		collide(other->getRightChild());
	} else {
		getLeftChild()->collide(other);
		getRightChild()->collide(other);
	}
}
