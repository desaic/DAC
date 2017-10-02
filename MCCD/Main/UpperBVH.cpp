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

#include "UpperBVH.h"
#include "DeformBVH.h"
#include "DeformModel.h"

static DeformModel *s_mdl;
static DeformBVHTree **s_trees;

UpperTree::UpperTree(DeformModel *mdl, bool ccd)
{
	_mdl = mdl;

	// need first decompose ...
	assert(mdl->_parts != NULL);

	_lowers = new DeformBVHTree *[mdl->_num_parts];
	for (unsigned int i=0; i<mdl->_num_parts; i++)
	{
		_lowers[i] = new DeformBVHTree(mdl, ccd, i);
		_lowers[i]->refit(false);
	}

	build_tree();
}

void
UpperTree::visulization(int level, bool upper)
{
	if (upper)
		getRoot()->visulization(level);
	else {
		for (unsigned int i=0; i<_mdl->_num_parts; i++)
			_lowers[i]->visulization(level, false);
	}
}

void
UpperTree::self_collide()
{
	s_trees = _lowers;
	s_mdl = _mdl;
	getRoot()->self_collide();
}

unsigned int g_changed;
extern unsigned int g_frame;

void
UpperTree::rebuild(bool ccd)
{
	g_changed = 0;

	if (_mdl->_num_parts >= 11) {
#pragma omp parallel for
		for (int i=0; i<_mdl->_num_parts; i++) {
			_lowers[i]->refit(false);
		}
	} else
		for (unsigned int i=0; i<_mdl->_num_parts; i++) {
			_lowers[i]->refit();
		}

	_aligned_free(_nodes);

	build_tree();

	s_trees = _lowers;
	getRoot()->refit();
}

void
UpperTree::refit(bool ccd)
{
	g_changed = 0;

	if (_mdl->_num_parts >= 11) {
#pragma omp parallel for
		for (int i=0; i<_mdl->_num_parts; i++) {
			_lowers[i]->refit(false);
		}
	} else
		for (unsigned int i=0; i<_mdl->_num_parts; i++) {
			_lowers[i]->refit();
			//_lowers[i]->check_convex();
		}

	s_trees = _lowers;
	getRoot()->refit();
}

//###################################################################

void
UpperNode::visulization(int level)
{
	if (isLeaf())
		_box.visulization();
	else
		if ((level > 0)) {
			if (level == 1)
				_box.visulization();
			else
			if (getLeftChild())
				getLeftChild()->visulization(level-1);
			if (getRightChild())
				getRightChild()->visulization(level-1);
		}
}

void
UpperNode::self_collide()
{
	if (isLeaf()) {
		s_trees[getID()]->self_collide();
		return;
	}

	getLeftChild()->self_collide();

	if (getRightChild()) {
		getRightChild()->self_collide();
		getLeftChild()->collide(getRightChild());
	}
}

void
UpperNode::collide(UpperNode *other)
{
	s_mdl->_counts[0]._num_box_tests++;
	if (!_box.overlaps(other->_box))
		return;

	if (isLeaf() && other->isLeaf()) {
		s_trees[getID()]->collide(s_trees[other->getID()]);
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

void
UpperNode::refit()
{
	if (isLeaf()) {
		_box = s_trees[getID()]->box();
	} else {
		getLeftChild()->refit();

		if (getRightChild())
		getRightChild()->refit();

		if (getRightChild())
			_box = getLeftChild()->_box + getRightChild()->_box;
		else
			_box = getLeftChild()->_box;
	}
}