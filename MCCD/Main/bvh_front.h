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

#include "box.h"
#include <vector>
#include <set>
using namespace std;

class bvh_front_list;
class bvh_front_node;
class non_adjacent_pair_list;

class DeformBVHNode;
class DeformModel;

typedef set<bvh_front_node> bvh_front_set;

class bvh_front_node {
	char _flag;
	DeformBVHNode *_left, *_right;

public:
	FORCEINLINE bvh_front_node(DeformBVHNode *left, DeformBVHNode *right)
	{
		if (left > right) {
			_left = left;
			_right = right;
		} else {
			_left = right;
			_right = left;
		}

		_flag = 0;
	}

	void update(bvh_front_list &, non_adjacent_pair_list &, int tid);

	bool operator == (const bvh_front_node &other) const {
		return	(_left == other._left && _right == other._right);
	}

	bool operator < (const bvh_front_node &other) const {
		if (_left == other._left)
			return _right < other._right;
		else
			return _left < other._left;
	}
};

class bvh_front_list : public vector<bvh_front_node> {
public:
	void propogate(DeformModel *mdl);
};
