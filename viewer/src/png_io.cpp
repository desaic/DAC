
#include "png_io.h"
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

#ifdef FREEIMAGE_AVAILABLE
#include "FreeImage.h"
void save_png(std::string & filename, BYTE * img, int w, int h)
{
  FIBITMAP* Image = FreeImage_ConvertFromRawBits(img, w, h, 3*w, 24, 0, 0, 0, false);
  FreeImage_Save(FIF_PNG, Image, filename.c_str(), 0);
  FreeImage_Unload(Image);
}
#else
void save_png(std::string & filename, BYTE * img, int w, int h)
{}
#endif
