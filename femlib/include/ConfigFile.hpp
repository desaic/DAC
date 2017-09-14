#ifndef CONFIGFILE_HPP
#define CONFIGFILE_HPP
#include <map>
#include <string>
#include <vector>

class Material;
class Quadrature;

class ConfigFile
{
public:
  ///@brief loads a config file.
  /// The absolute path of the directory is also stored in dir.
  int load(const std::string filename);

  void getFloat(const std::string & key, float & val)const;

  ///@return true of val is written into. Otherwise val is untouched.
  bool getInt(const std::string & key, int &val )const;
  bool getIntVector(const std::string & key, std::vector<int> &a )const;

  ///@param val not modified if key not present.
  void getBool(const std::string & key, bool & val)const;
  std::vector<float> getFloatVector(const std::string&key)const;
  ///@return empty string if key is not present.
  std::string getString(const std::string & key)const;
  std::vector<std::string> getStringVector(const std::string & key)const;

  bool hasOpt(const std::string & key)const;

  ///@brief name of directory containing the config file
  std::string dir;
private:
  std::map<std::string, std::vector<std::string> > m;
};

/// \brief loads a list of materials from "materialfile" option in config file.
int loadMaterials(std::vector<Material*> & materials, const ConfigFile & conf, const Quadrature* quadrature);
#endif
