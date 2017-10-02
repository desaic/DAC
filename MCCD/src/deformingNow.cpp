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

#include "DeformModel.h"
#include <stdio.h>
#include "timing.h"

CBVHTimer tm;

#define frand48()  ((((float) rand()) / ((float) RAND_MAX)))

static DeformModel *mdl;

void initModel(double * verts0, double * verts1, int nv, int * trigs, int nt, int num_frame, bool ccd)
{
  TIMING_BEGIN
    mdl = new DeformModel(verts0, nv, trigs, nt, 1.0f);
  TIMING_END("Load models")

  TIMING_BEGIN
    mdl->Decompose();
  TIMING_END("Decompose")

  TIMING_BEGIN
    mdl->BuildBVH(ccd);
  TIMING_END("Build BVH")
}

void quitModel()
{
	delete mdl;
}

void drawModel(int level, bool upper)
{
	mdl->Display();
	mdl->DisplayBVH(level, upper);

}

void dynamicModel_AABB(int t, int circle, bool refit, bool ccd)
{
	tm.startTiming(9);
	tm.startTiming(10);

	mdl->RebuildBVH(ccd);

	tm.endTiming(10);
	tm.endTiming(9);

	mdl->ResetCounter();

	tm.startTiming(0);

	mdl->SelfCollide(ccd);
	
  //mdl->UpdateCollide();

	tm.endTiming(0);

	tm.incRecord(mdl->NumBoxTest(), mdl->NumVFTest(), mdl->NumEETest(), mdl->NumVFTrue(), mdl->NumEETrue());
}

extern bool *pb;
static double g_total = 0;
//extern void endCapture();

void dynamicModel(bool refit, bool ccd, bool obb)
{
	double ttmp = omp_get_wtime();
  
  mdl->Deform();
	
  tm.updatTiming();

	dynamicModel_AABB(0, 1, refit, ccd);

	g_total += omp_get_wtime() - ttmp;
  printf("total: %g seconds\n", g_total);
  tm.report();
  //endCapture();
}
