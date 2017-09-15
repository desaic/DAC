#ifndef SAMPLEFORCES_HPP
#define SAMPLEFORCES_HPP

#include <vector>
#include "ConfigFile.hpp"

class Material;

void SampleForces(const ConfigFile & conf);
void fitParams(const ConfigFile & conf);
void extractComb(const ConfigFile & conf);
void makeCoarseMats(std::vector<std::vector<double > > & params,
  std::vector<Material*> & materials);
#endif