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

#include "logger.h"
#include "TrigMesh.hpp"
#include <fstream>
Logger *logger;

char DATA_PATH[512];
int NUM_FRAME;

float MODEL_SCALE = 1.f;
float DISP_SCALE = 1.f;
int SLICES = 5;
int START_FRAME=0;//320;
bool CCD = true;

static int s_count = 0;
extern void dynamicModel(int t, int start, int circle, bool refit, bool ccd);
extern void initModel(char *, int, float, bool ccd);

void initModel(double * verts, int nv, int * trigs, int nt, int num_frame, float scale, bool ccd);

void endCapture()
{
	NULL;
}

int main(int argc, char **argv)
{
	if (argc<2) {
		printf("Usage: %s obj1 obj2\n", argv[0]);
		exit(0);
	}
  TrigMesh tm[2];
  for (int i = 0; i < 2; i++) {
    std::ifstream in(argv[i+1]);
    if (!in.good()) {
      std::cout << "Can't open " << argv[i + 1] << "\n";
      return -1;
    }
    tm[i].load(in);
    in.close();
  }
  int NUM_FRAME = 2;
  std::cout << "size of a v and t" << sizeof(tm[0].v[0]) << " " << sizeof(tm[0].t[0]) << "\n";
	initModel((double*)tm[0].v.data(), tm[0].v.size(), (int*)tm[0].t.data(), tm[0].t.size(), 
    2, MODEL_SCALE, CCD);
  //initModel("", 2, 1.0, true);
	dynamicModel(0, 0, 1, true, CCD);
	
	return 0;
}
