#ifndef PNG_IO_H
#define PNG_IO_H

#include <string>
typedef unsigned char BYTE;
//not thread safe. Can be modified to be if needed.
void save_png(std::string & filename, BYTE * img, int w, int h);

#endif