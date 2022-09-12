#include "OutputBuffer.hpp"
#include <cassert>
#include <cstring> // memcpy

void outputBuffer::clear()
{
  this->states.clear();
  currentReadIter = 0;
}

void outputBuffer::writeNewStep()
{
  this->states.push_back(State());
}

bool outputBuffer::canRead()
{
  return currentReadIter < this->states.size();
}
void outputBuffer::readNextStep()
{
  ++currentReadIter;
}

double *outputBuffer::getDoubleData(const std::string &name)
{
  assert(states.size() > currentReadIter);
  auto &data_vector = this->states[currentReadIter].data_double[name];
  return data_vector.data();
}

ITG *outputBuffer::getITGData(const std::string &name)
{
  assert(states.size() > currentReadIter);
  auto &data_vector = this->states[currentReadIter].data_itg[name];
  return data_vector.data();
}

void outputBuffer::writeDoubleData(const std::string &name, double *data, unsigned n)
{
  assert(!states.empty());
  this->states.back().data_double[name] = std::vector<double>(data, data + n);
}

void outputBuffer::writeITGData(const std::string &name, ITG *data, unsigned n)
{
  this->states.back().data_itg[name] = std::vector<ITG>(data, data + n);
}

/* Implement the C interface */
outputBuffer *BufferCreate()
{
  return new outputBuffer();
}
void BufferFree(outputBuffer *buffer)
{
  delete buffer;
}
void BufferSaveDouble(outputBuffer *buffer, const char *name, double *data, unsigned length)
{
  buffer->writeDoubleData(name, data, length);
}
void BufferSaveITG(outputBuffer *buffer, const char *name, ITG *data, unsigned length)
{
  buffer->writeITGData(name, data, length);
}
void BufferLoadDouble(outputBuffer *buffer, const char *name, double *data, unsigned length)
{
  memcpy(data, buffer->getDoubleData(name), length); //Length from the buffer or from ptr ?
}
void BufferLoadITG(outputBuffer *buffer, const char *name, ITG *data, unsigned length)
{
  memcpy(data, buffer->getITGData(name), length); //Length from the buffer or from ptr ?
}
int BufferCanRead(outputBuffer *buffer)
{
  return buffer->canRead();
}
void BufferReadNext(outputBuffer *buffer)
{
  buffer->readNextStep();
}
void BufferWriteNewStep(outputBuffer *buffer)
{
  buffer->writeNewStep();
}
void BufferClear(outputBuffer *buffer)
{
  buffer->clear();
}
unsigned BufferStoredStates(outputBuffer *buffer)
{
  return buffer->states.size();
}