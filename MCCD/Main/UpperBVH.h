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

class DeformModel;
class DeformBVHTree;

class UpperNode {
	BOX  _box;

	int _child;

public:
	UpperNode();

	void Construct(unsigned int);
	void Construct(unsigned int *, unsigned int);

	~UpperNode();

	void collide(UpperNode *);
	void self_collide();
	void visulization(int level);
	void refit();

	FORCEINLINE UpperNode *getLeftChild() { return this - _child; }
	FORCEINLINE UpperNode *getRightChild() { return this - _child + 1; }

	FORCEINLINE int getID() { return _child; }
	FORCEINLINE int isLeaf() { return _child >= 0; }

friend class UpperTree;
};

class UpperTree {
	UpperNode *_nodes;

	DeformModel	*_mdl;
	DeformBVHTree **_lowers;

	void build_tree();

public:
	UpperTree(DeformModel *, bool);
	~UpperTree();

	void Construct(DeformModel *);
	void ConstructCCD(DeformModel *);

	void visulization(int level, bool upper);
	void self_collide();

	void refit(bool ccd);
	void rebuild(bool ccd);

	FORCEINLINE UpperNode *getRoot() { return _nodes; }
};
