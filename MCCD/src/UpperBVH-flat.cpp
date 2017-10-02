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
#include "aap.h"
static DeformBVHTree **s_trees;

#include "timing.h"
extern CBVHTimer tm;

static UpperNode *s_current_node;

void
UpperTree::build_tree()
{
	int count = _mdl->_num_parts;

	// now we build _root
	BOX total;
	for (int i=0; i<count; i++)
		total += _lowers[i]->box();

	aap pln(total);
	unsigned int *idx_buffer = new unsigned int[count];
	unsigned int left_idx = 0, right_idx = count;

	for (int i=0; i<count; i++) {
		if (pln.inside(_lowers[i]->box().center()))
			idx_buffer[left_idx++] = i;
		else
			idx_buffer[--right_idx] = i;
	}

	//_root = new UpperNode();
	//_root->_box = total;
	//_nodes = new UpperNode[count*2-1];
	_nodes = (UpperNode*) _aligned_malloc((count*2-1) * sizeof(UpperNode), 16);
	_nodes[0]._box = total;
	s_current_node = _nodes+3;

	s_trees = _lowers;

	int hal = count/2;
	if (count == 1) {
		_nodes[0]._child = 0;
	} else
	if (left_idx == 0 || left_idx == count)
	{
		_nodes[0]._child = -1;
		_nodes[0].getLeftChild()->Construct(idx_buffer, hal);
		_nodes[0].getRightChild()->Construct(idx_buffer+hal, count-hal);
	}
	else {
		_nodes[0]._child = -1;
		_nodes[0].getLeftChild()->Construct(idx_buffer, left_idx);
		_nodes[0].getRightChild()->Construct(idx_buffer+left_idx, count-left_idx);
	}

	delete [] idx_buffer;
}

UpperTree::~UpperTree()
{
	//delete [] _nodes;
	_aligned_free(_nodes);

	for (unsigned int i=0; i<_mdl->_num_parts; i++)
		delete _lowers[i];

	delete [] _lowers;
}

void
UpperNode::Construct(unsigned int id)
{
	_child = id;
}

void
UpperNode::Construct(unsigned int *lst, unsigned int lst_num)
{
	if (lst_num == 1){
		_child = lst[0];
		_box = s_trees[getID()]->box();
		return;
	}

	_box.empty();
	//try to split them
	for (unsigned int t=0; t<lst_num; t++) {
		int i=lst[t];
		_box += s_trees[i]->box();
	}

	if (lst_num == 2) {// must split it!
		getLeftChild()->Construct(lst[0]);
		getRightChild()->Construct(lst[1]);
		return;
	}

	aap pln(_box);
	unsigned int left_idx=0, right_idx=lst_num-1;
	
	for (unsigned int t=0; t<lst_num; t++) {
		int i=lst[left_idx];
		if (pln.inside(s_trees[i]->box().center()))
			left_idx++;
		else {// swap it
			unsigned int tmp = i;
			lst[left_idx] = lst[right_idx];
			lst[right_idx--] = tmp;
		}
	}

	int hal = lst_num/2;
	if (left_idx == 0 || left_idx == lst_num)
	{
		getLeftChild()->Construct(lst, hal);
		getRightChild()->Construct(lst+hal, lst_num-hal);
	}
	else {
		getLeftChild()->Construct(lst, left_idx);
		getRightChild()->Construct(lst+left_idx, lst_num-left_idx);
	}
}


UpperNode::UpperNode()
{
	_child = 0;
}

UpperNode::~UpperNode()
{
	NULL;
}
