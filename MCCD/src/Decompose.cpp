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

#include <windows.h>

#include <math.h>
#include <stdio.h>
#include <assert.h>
#include "DeformModel.h"

void
DeformModel::Decompose()
{
	_num_parts = 0;
	_parts = new unsigned int[_num_tri];
	for (unsigned int i=0; i<_num_tri; i++)
		_parts[i] = -1;

	bool *flags;
	flags = new bool[_num_vtx];

	bool done = false;
	unsigned int cur_part = 0;

	while (!done) {
		done = true;

		for (unsigned int i=0; i<_num_vtx; i++)
			flags[i] = false;

		bool first = true;
		bool find = true;

		while (find) {
			find = false;

			for (unsigned int i=0; i<_num_tri; i++) {
				if (first && _parts[i] == -1) {
					first = false;
					done = false;
					find = true;

					_parts[i] = cur_part;

					tri3f *tri = _tris+i;
					for (int j=0; j<3; j++)
						flags[tri->id(j)] = true;
				}

				if (_parts[i] == -1) {
					tri3f *tri = _tris+i;

					bool b=false;
					for (int j=0; j<3; j++)
						if (flags[tri->id(j)] == true) {
							find = true;
							b = true;
							break;
						}

					if (b) {
						_parts[i] = cur_part;
						for (int j=0; j<3; j++)
							flags[tri->id(j)] = true;
					}
				}
			}
		}

		cur_part++;
	}

	_num_parts = cur_part-1;

	if (_colors == NULL) {
		static unsigned char iColor[] = {
			255, 128, 0,
			0, 255, 0,
			0, 0, 255,
			255, 0, 0,
			255, 255, 0,
			0, 255, 255,
			255, 0, 255,
			128, 128, 0,
			0, 255, 128,
			128, 0, 0,
			128, 255, 0,
			255, 128, 0,
			0, 128, 255,
			0, 255, 0,
			0, 255, 0,
			0, 255, 0,
			0, 255, 0,
		};
		static int length = sizeof(iColor)/3;

		_colors = new color3[_num_vtx];

		for (unsigned int j=0; j<_num_tri; j++)
		{
			for (int i=0; i<3; i++) {
				unsigned int vtx_id = _tris[j].id(i);
				_colors[vtx_id].set(iColor+(_parts[j]%length)*3);
			}
		}
	}
}