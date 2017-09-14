#ifndef FILEUTIL_HPP
#define FILEUTIL_HPP
#include <vector>
#include <iostream>
#include <fstream>
#include "Eigen/Dense"
class FileUtilIn
{
public:
  FileUtilIn();
  FileUtilIn(const std::string filename);
  void open(const std::string filename, std::ios_base::openmode mode= std::ios_base::in);
  void close();
  void readFArr(std::vector<float> & arr);
  void readFArr(std::vector<float> & arr, int size);
  void readArr(std::vector<Eigen::Vector3f> & arr);
  bool good(){ return in.good(); }
  template <typename T>
  friend std::istream & operator>>(std::istream & ii, std::vector<T> & arr);
  std::ifstream in;
};

class FileUtilOut
{
public:

  FileUtilOut();
  FileUtilOut(const std::string filename);
  void open(const std::string filename, std::ios_base::openmode mode = std::ios_base::out);
  void close();

  template <typename T>
  friend std::ostream & operator<<(std::ostream & o, const std::vector<T> & arr);

  void writeArr(const std::vector<Eigen::Vector3f> & arr);
  std::ofstream out;
};

///@brief make a string dirname/prefixseq.txt
std::string
sequenceFilename(int seq, const char * dirname = 0,
                          const char * prefix  = 0,
                 int padding = 0, const char * suffix= 0);

std::string directoryName(std::string filename);

template <typename T>
std::vector<std::vector<T> > loadArr2d(std::string filename)
{
  FileUtilIn in(filename);
  std::vector<std::vector<T> > arr;
  if (!in.good()){
    return arr;
  }
  int nrow, ncol;
  in.in >> nrow >> ncol;
  arr.resize(nrow);
  for (int i = 0; i < nrow; i++){
    arr[i].resize(ncol);
    for (int j = 0; j < ncol; j++){
      in.in >> arr[i][j];
    }
  }
  return arr;
}
#endif // FILEUTIL_HPP
