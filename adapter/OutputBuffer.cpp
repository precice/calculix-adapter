/**********************************************************************************************
 *                                                                                            *
 *       CalculiX adapter for heat transfer coupling and mechanical FSI using preCICE         *
 *       Heat transfer adapter developed by Lucía Cheung with the support of SimScale GmbH    *
 *                                                                                            *
 *       Adapter extended to fluid-structure interaction by Alexander Rusch                   *
 *                                                                                            *
 *********************************************************************************************/

#include "OutputBuffer.hpp"

void Buffer::clear() {
    this->states.clear();
    currentReadIter = 0;
}

void Buffer::writeNewStep() {
    this->states.push_back(State());
}

void Buffer::startRead(){
    currentReadIter = 0;
}

bool Buffer::readNext() {
    ++currentReadIter;
    return currentReadIter < this->states.size();
}

double* Buffer::getDoubleData(const std::string &name) {
    return this->states[currentReadIter].data_double[name].data();
}

ITG* Buffer::getITGData(const std::string &name) {
    return this->states[currentReadIter].data_itg[name].data();
}

unsigned Buffer::getDoubleDataSize(const std::string &name) {
    return this->states[currentReadIter].data_double[name].size();
}

unsigned Buffer::getITGDataSize(const std::string &name) {
    return this->states[currentReadIter].data_itg[name].size();
}

void Buffer::writeDoubleData(const std::string &name, double *data, unsigned n) {
    this->states.back().data_double[name] = std::vector<double> (data, data + n);
}

void Buffer::writeITGData(const std::string &name, ITG *data, unsigned n) {
    this->states.back().data_itg[name] = std::vector<ITG> (data, data + n);
}