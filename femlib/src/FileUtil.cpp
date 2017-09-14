#include "FileUtil.hpp"
#include <stdlib.h>
#include <string>
#include <iomanip>
#include <sstream>

std::istream & operator>> (std::istream  & in, Eigen::Vector3f & v)
{
  for(int ii= 0; ii<3; ii++){
    in>>v[ii];
  }
  return in;
}

void
FileUtilIn::readArr(std::vector<Eigen::Vector3f> & arr)
{
  for(unsigned int ii = 0;ii<arr.size();ii++){
    in >> arr[ii];
  }
}

void
FileUtilIn::open(const std::string filename,std::ios_base::openmode mode )
{
  in.open(filename,mode);
  if(!in.good()){
    std::cout<<"cannot open "<<filename<<"\n";
  }
}

void FileUtilIn::close()
{
  in.close();
}

void
FileUtilIn::readFArr(std::vector<float> & arr)
{
  int size;
  in>>size;
  arr.resize(size);
  for(int ii = 0;ii<size;ii++){
    in>>arr[ii];
  }
}

void
FileUtilIn::readFArr(std::vector<float> & arr, int size)
{
  arr.resize(size);
  for(int ii = 0;ii<size;ii++){
    in>>arr[ii];
  }
}

template <typename T>
std::istream &
operator>>(std::istream & i, std::vector<T> & arr)
{
  int size = 0;
  i>>size;
  arr.resize(size);
  for(int ii = 0;ii<size;ii++){
    i>>arr[ii];
  }
}

FileUtilIn::FileUtilIn(const std::string filename):in(filename)
{
  if(!in.good()){
    std::cout<<"cannot open "<<filename<<"\n";
  }
}

FileUtilIn::FileUtilIn()
{}

void
FileUtilOut::writeArr(const std::vector<Eigen::Vector3f> & arr)
{
  for(unsigned int ii = 0;ii<arr.size();ii++){
    out<<arr[ii]<<"\n";
  }
}

void
FileUtilOut::open(const std::string filename, std::ios_base::openmode mode)
{
  out.open(filename,mode);
  if(!out.good()){
    std::cout<<"cannot open output"<<filename<<"\n";
  }
}

void
FileUtilOut::close()
{
  out.close();
}

FileUtilOut::FileUtilOut(const std::string filename):out(filename)
{
  if(!out.good()){
    std::cout<<"cannot open output"<<filename<<"\n";
  }
}

FileUtilOut::FileUtilOut()
{}

template <typename T>
std::ostream &
operator<<(std::ostream & o, const std::vector<T> & arr)
{
  o<<arr.size()<<"\n";
  for(unsigned int ii = 0 ;ii<arr.size();ii++){
    o<<arr[ii]<<"\n";
  }
}

std::string
sequenceFilename(int seq, const char * dirname,
                          const char * prefix,
                 int padding, const char * _suffix)
{
  std::stringstream ss;
  if(dirname != 0){
    ss<<dirname<<"/";
  }
  if(prefix != 0){
    ss<<prefix;
  }
  std::string suffix(".txt");
  if(_suffix!=0){
    suffix = std::string(_suffix);
  }
  if(padding>0){
    ss << std::setfill('0') << std::setw(padding)<<seq<<suffix;
  }else{
    ss << seq<<suffix;
  }
  return ss.str();
}

std::string directoryName(std::string filename)
{
  char * fullpath;
  std::string path;
  std::string s;

#ifdef __linux__
  char resolve[PATH_MAX];
  fullpath = realpath(filename.c_str(), resolve);
  if (fullpath == NULL){
    std::cout << "Error resolving " << filename << "\n";
  }
  s = std::string(fullpath);
#elif _WIN32
  char full[_MAX_PATH];
  if (_fullpath(full, filename.c_str(), _MAX_PATH) != NULL){
    s = std::string(full);
  }
#endif
  std::size_t found = s.find_last_of("/\\");
  path = s.substr(0, found);
  path += "/";
  return path;
}
