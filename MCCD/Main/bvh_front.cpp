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

#include "bvh_front.h"
#include "DeformBVH.h"
#include "DeformModel.h"
#include <omp.h>

static DeformModel *s_mdl1 = NULL;

#include "tri_pair.h"
extern non_adjacent_pair_list non_adj_list;

#include "timing.h"
extern CBVHTimer tm;

void
bvh_front_node::update(bvh_front_list &appended, non_adjacent_pair_list &result, int tid)
{
	if (_flag != 0) return;

	if (_left->isLeaf() && _right->isLeaf()) {
		s_mdl1->_counts[tid]._num_box_tests++;
		if (_left->_box.overlaps(_right->_box)) {
			// still intersecting, so keep it in front, and record related intersection
			result.push_back(non_adjacent_pair(_left->getTriID(), _right->getTriID()));

		}
	} else {
		s_mdl1->_counts[tid]._num_box_tests++;
		if (_left->_box.overlaps(_right->_box)) { // need to expand here
				_flag = 1; // this node is invaild

				if (_left->isLeaf()) {
					_left->append_collide(_right->getLeftChild(), appended, result, tid);
					_left->append_collide(_right->getRightChild(), appended, result, tid);
				} else {
					_left->getLeftChild()->append_collide(_right, appended, result, tid);
					_left->getRightChild()->append_collide(_right, appended, result, tid);
				}
		}
	}
}

void
bvh_front_list::propogate(DeformModel *mdl)
{
	s_mdl1 = mdl;

	bvh_front_list append[16];
	non_adjacent_pair_list result[16];

	int length = size();
	int frag_length = length/16;

#pragma omp parallel sections
	{
 #pragma omp section
			{
				for (int i=0; i<frag_length; i++)
					(*this)[i].update(append[0], result[0], 0);
			}
 #pragma omp section
			{
				for (int i=frag_length; i<frag_length*2; i++)
					(*this)[i].update(append[1], result[1], 1);
			}
 #pragma omp section
			{
				for (int i=frag_length*2; i<frag_length*3; i++)
					(*this)[i].update(append[2], result[2], 2);
			}
 #pragma omp section
			{
				for (int i=frag_length*3; i<frag_length*4; i++)
					(*this)[i].update(append[3], result[3], 3);
			}
 #pragma omp section
			{
				for (int i=frag_length*4; i<frag_length*5; i++)
					(*this)[i].update(append[4], result[4], 4);
			}
 #pragma omp section
			{
				for (int i=frag_length*5; i<frag_length*6; i++)
					(*this)[i].update(append[5], result[5], 5);
			}
 #pragma omp section
			{
				for (int i=frag_length*6; i<frag_length*7; i++)
					(*this)[i].update(append[6], result[6], 6);
			}
 #pragma omp section
			{
				for (int i=frag_length*7; i<frag_length*8; i++)
					(*this)[i].update(append[7], result[7], 7);
			}
 #pragma omp section
			{
				for (int i=frag_length*8; i<frag_length*9; i++)
					(*this)[i].update(append[8], result[8], 8);
			}
 #pragma omp section
			{
				for (int i=frag_length*9; i<frag_length*10; i++)
					(*this)[i].update(append[9], result[9], 9);
			}
 #pragma omp section
			{
				for (int i=frag_length*10; i<frag_length*11; i++)
					(*this)[i].update(append[10], result[10], 10);
			}
 #pragma omp section
			{
				for (int i=frag_length*11; i<frag_length*12; i++)
					(*this)[i].update(append[11], result[11], 11);
			}
 #pragma omp section
			{
				for (int i=frag_length*12; i<frag_length*13; i++)
					(*this)[i].update(append[12], result[12], 12);
			}
 #pragma omp section
			{
				for (int i=frag_length*13; i<frag_length*14; i++)
					(*this)[i].update(append[13], result[13], 13);
			}
 #pragma omp section
			{
				for (int i=frag_length*14; i<frag_length*15; i++)
					(*this)[i].update(append[14], result[14], 14);
			}
 #pragma omp section
			{
				for (int i=frag_length*15; i<length; i++)
					(*this)[i].update(append[15], result[15], 15);
			}
 	}

	for (int i=0; i<16; i++) {
		insert(end(), append[i].begin(), append[i].end());
		non_adj_list.insert(non_adj_list.end(), 	result[i].begin(), result[i].end());
	}
}


