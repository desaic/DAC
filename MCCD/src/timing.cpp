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

# include "timing.h"
# include <stdio.h>
# include <assert.h>

void
CBVHTimer::resetTiming()
{
	_count = 0;

	for (int i=0; i<10; i++) {
		_rets[i] =  0;
	}

	_box_test = _vf_true = _ee_true = _vf_test = _ee_test = 0L;

}

void
CBVHTimer::incRecord(int box_test, int vf_test, int ee_test, int vf_true, int ee_true)
{
	_box_test += box_test;
	_vf_test += vf_test;
	_ee_test += ee_test;
	_vf_true += vf_true;
	_ee_true += ee_true;
}

void
CBVHTimer::report()
{
	printf("frame %d: =================================\n", _count);
	printf("intersection time: %g\n", _rets[0]);


	printf("    updating front & collecting: %g\n", _rets[12]);
	printf("    collecting nonadjacent: %g\n", _rets[1]);
	printf("    processing nonadjacent: %g\n", _rets[4]);

	printf("    processing adjacent: %g\n", _rets[2]);
	printf("    merging bvh_front_list: %g\n", _rets[5]);

	printf("BVH and update time: %g\n", _rets[9]);
	printf("    refit BVH: %g\n", _rets[10]);

	printf("    update vtx, boxes for V/E/F: %g\n", _rets[11]);
	printf("all vf tests: %ld, true vf tests: %ld, ratio: %f\n", _vf_test, _vf_true, float(_vf_true)/_vf_test);
	printf("all ee tests: %ld, true ee tests: %ld, ratio: %f\n", _ee_test, _ee_true, float(_ee_true)/_ee_test);
	printf("box test: %ld, contact: %ld\n", _box_test, _contact);
}
