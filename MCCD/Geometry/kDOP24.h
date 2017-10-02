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

#include "vec3f.h"
#include <float.h>
#include <stdlib.h>

#define MAX(a,b)	((a) > (b) ? (a) : (b))
#define MIN(a,b)	((a) < (b) ? (a) : (b))

inline void MinMax(float p, float &mi, float &ma)
{
	if (p > ma) ma = p;
	if (p < mi) mi = p;
}

inline void MinMax(float a, float b, float &mi, float &ma)
{
	if (a > b)
		mi = b, ma = a;
	else
		mi = a, ma = b;
}

/*
     * k=18 (Aabb + 12 diagonal planes that "cut off" some space of the edges):
     * \li (-1,0,0) and (1,0,0)  -> indices 0 and 9
     * \li (0,-1,0) and (0,1,0)  -> indices 1 and 10
     * \li (0,0,-1) and (0,0,1)  -> indices 2 and 11
     * \li (-1,-1,0) and (1,1,0) -> indices 3 and 12
     * \li (-1,0,-1) and (1,0,1) -> indices 4 and 13
     * \li (0,-1,-1) and (0,1,1) -> indices 5 and 14
     * \li (-1,1,0) and (1,-1,0) -> indices 6 and 15
     * \li (-1,0,1) and (1,0,-1) -> indices 7 and 16
     * \li (0,-1,1) and (0,1,-1) -> indices 8 and 17
*/

class kDOP24 {
public:
	FORCEINLINE static void getDistances(const vec3f& p, float d[])
	{
		d[0] = p[0] + p[1];
		d[1] = p[0] + p[2];
		d[2] = p[1] + p[2];
		d[3] = p[0] - p[1];
		d[4] = p[0] - p[2];
		d[5] = p[1] - p[2];

		d[6] = p[0]+p[1]-p[2];
		d[7] = p[0]+p[2]-p[1];
		d[8] = p[1]+p[2]-p[0];
	}

public:
	float _dist[24];

	FORCEINLINE kDOP24() {
		empty();
	}

	FORCEINLINE kDOP24(const vec3f &v) {
		for (int i=0; i<3; i++)
			_dist[i] = _dist[12+i] = v[i];

		float d[9];
		getDistances(v, d);
		for (int i=0; i<9; i++)
			_dist[3+i] = _dist[15+i] = d[i];
 	}

	FORCEINLINE kDOP24(const vec3f &a, const vec3f &b) {
		for (int i=0; i<3; i++)
			MinMax(a[i], b[i], _dist[i], _dist[12+i]);

		float ad[9], bd[9];
		getDistances(a, ad);
		getDistances(b, bd);
		for (int i=0; i<9; i++)
			MinMax(ad[i], bd[i], _dist[3+i], _dist[15+i]);
	}

	FORCEINLINE bool overlaps(const kDOP24& b) const
	{
		for (int i=0; i<12; i++) {
			if (_dist[i] > b._dist[i+12]) return false;
			if (_dist[i+12] < b._dist[i]) return false;
		}

		return true;
	}

	FORCEINLINE kDOP24 &operator += (const vec3f &p)
	{
		for (int i=0; i<3; i++)
			MinMax(p[i], _dist[i], _dist[12+i]);

		float pd[9];
		getDistances(p, pd);
		for (int i=0; i<9; i++)
			MinMax(pd[i], _dist[3+i], _dist[15+i]);

		return *this;
	}

	FORCEINLINE kDOP24 &operator += (const kDOP24 &b)
	{
		for (int i=0; i<12; i++) {
			_dist[i]  = MIN(b._dist[i], _dist[i]);
			_dist[i+12]  = MAX(b._dist[i+12], _dist[i+12]);
		}
		return *this;
	}

	FORCEINLINE kDOP24 operator + ( const kDOP24 &v) const
	{ kDOP24 rt(*this); return rt += v; }

	FORCEINLINE float width()  const { return _dist[12] - _dist[0]; }
	FORCEINLINE float height() const { return _dist[13] - _dist[1]; }
	FORCEINLINE float depth()  const { return _dist[14] - _dist[2]; }
	FORCEINLINE float volume() const { return width()*height()*depth(); }

	FORCEINLINE vec3f center() const { 
		return vec3f(_dist[0]+_dist[12], _dist[1]+_dist[13], _dist[2]+_dist[14])*0.5f;
	}

	FORCEINLINE void empty() {
		for (int i=0; i<12; i++) {
			_dist[i] = FLT_MAX;
			_dist[i+12] = -FLT_MAX;
		}
	}

	void visulization() {}
};