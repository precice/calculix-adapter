/**********************************************************************************************
 *                                                                                            *
 *       CalculiX adapter for heat transfer coupling and mechanical FSI using preCICE         *
 *       Heat transfer adapter developed by Lucía Cheung with the support of SimScale GmbH    *
 *                                                                                            *
 *       Adapter extended to fluid-structure interaction by Alexander Rusch                   *
 *                                                                                            *
 *********************************************************************************************/

#include "OutputBuffer.hpp"
#include <cstring> // memcpy

void outputBuffer::clear() {
    this->states.clear();
    currentReadIter = 0;
}

void outputBuffer::writeNewStep() {
    this->states.push_back(State());
}

void outputBuffer::startRead(){
    currentReadIter = 0;
}

bool outputBuffer::readNext() {
    ++currentReadIter;
    return currentReadIter < this->states.size();
}

double* outputBuffer::getDoubleData(const std::string &name) {
    return this->states[currentReadIter].data_double[name].data();
}

ITG* outputBuffer::getITGData(const std::string &name) {
    return this->states[currentReadIter].data_itg[name].data();
}

unsigned outputBuffer::getDoubleDataSize(const std::string &name) {
    return this->states[currentReadIter].data_double[name].size();
}

unsigned outputBuffer::getITGDataSize(const std::string &name) {
    return this->states[currentReadIter].data_itg[name].size();
}

void outputBuffer::writeDoubleData(const std::string &name, double *data, unsigned n) {
    this->states.back().data_double[name] = std::vector<double> (data, data + n);
}

void outputBuffer::writeITGData(const std::string &name, ITG *data, unsigned n) {
    this->states.back().data_itg[name] = std::vector<ITG> (data, data + n);
}

/* Implement the C interface */
outputBuffer* BufferCreate() {
    return new outputBuffer();
}
void BufferFree(outputBuffer* buffer) {
    delete buffer;    
}
void BufferSaveDouble(outputBuffer* buffer, const char * name, double * data, unsigned length) {
    buffer->writeDoubleData(name, data, length);
}
void BufferSaveITG(outputBuffer* buffer, const char * name, ITG * data, unsigned length) {
    buffer->writeITGData(name, data, length);
}
unsigned BufferGetLengthDouble(outputBuffer* buffer, const char * name) {
    return buffer->getDoubleDataSize(name);
}
unsigned BufferGetLengthITG(outputBuffer* buffer, const char * name) {
    return buffer->getITGDataSize(name);
}
void BufferLoadDouble(outputBuffer* buffer, const char * name, double * data, unsigned length) {
    memcpy(data, buffer->getDoubleData(name), length); //Length from the buffer or from ptr ?
}
void BufferLoadITG(outputBuffer* buffer, const char * name, ITG * data, unsigned length) {
    memcpy(data, buffer->getITGData(name), length); //Length from the buffer or from ptr ?
}
void BufferNextIter(outputBuffer* buffer) {
    buffer->readNext();
}
void BufferClear(outputBuffer* buffer) {
    buffer->clear();
}