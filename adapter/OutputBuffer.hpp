/*
  Implementation of the API defined in OutputBuffer.h
*/

#ifndef CCX_OUTPUT_BUFFER_HPP
#define CCX_OUTPUT_BUFFER_HPP

extern "C" {

#include "OutputBuffer.h"
}

#include <map>
#include <string>
#include <vector>

struct State {
  std::map<std::string, std::vector<double>> data_double;
  std::map<std::string, std::vector<ITG>>    data_itg;
};

class outputBuffer {
public:
  std::vector<State> states;

  void clear();
  void writeNewStep();
  void readNext();
  bool canRead();

  std::vector<double> &loadDouble(const std::string &name);
  std::vector<ITG> &   loadITG(const std::string &name);

  void saveDouble(const std::string &name, double *data, unsigned n);
  void saveITG(const std::string &name, ITG *data, unsigned n);

private:
  unsigned currentReadIter{0};
};

#endif
