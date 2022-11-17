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
void outputBuffer::readNext()
{
  ++currentReadIter;
}

std::vector<double> &outputBuffer::loadDouble(const std::string &name)
{
  assert(states.size() > currentReadIter);
  return this->states[currentReadIter].data_double[name];
}

std::vector<ITG> &outputBuffer::loadITG(const std::string &name)
{
  assert(states.size() > currentReadIter);
  return this->states[currentReadIter].data_itg[name];
}

void outputBuffer::saveDouble(const std::string &name, double *data, unsigned n)
{
  assert(!states.empty());
  this->states.back().data_double[name] = std::vector<double>(data, data + n);
}

void outputBuffer::saveITG(const std::string &name, ITG *data, unsigned n)
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
  buffer->saveDouble(name, data, length);
}
void BufferSaveITG(outputBuffer *buffer, const char *name, ITG *data, unsigned length)
{
  buffer->saveITG(name, data, length);
}
void BufferLoadDouble(outputBuffer *buffer, const char *name, double *data)
{
  const auto &vector = buffer->loadDouble(name);
  memcpy(data, vector.data(), vector.size() * sizeof(double));
}
void BufferLoadITG(outputBuffer *buffer, const char *name, ITG *data)
{
  const auto &vector = buffer->loadITG(name);
  memcpy(data, vector.data(), vector.size() * sizeof(ITG));
}
int BufferCanRead(outputBuffer *buffer)
{
  return buffer->canRead();
}
void BufferReadNext(outputBuffer *buffer)
{
  buffer->readNext();
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