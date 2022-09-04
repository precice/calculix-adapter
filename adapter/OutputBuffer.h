/**********************************************************************************************
 *                                                                                            *
 *       CalculiX adapter for heat transfer coupling and mechanical FSI using preCICE         *
 *       Heat transfer adapter developed by Lucía Cheung with the support of SimScale GmbH    *
 *                                                                                            *
 *       Adapter extended to fluid-structure interaction by Alexander Rusch                   *
 *                                                                                            *
 *********************************************************************************************/

#ifndef CCX_OUTPUTBUFFER_H
#define CCX_OUTPUTBUFFER_H

// To avoid including CalculiX.h
#ifdef LONGLONG
#define ITG long long
#define ITGFORMAT "lld"
#else
#define ITG int
#define ITGFORMAT "d"
#endif


struct outputBuffer;

outputBuffer* BufferCreate();
void BufferFree(outputBuffer* buffer);
void BufferSaveDouble(outputBuffer* buffer, const char * name, double * data, unsigned length);
void BufferSaveITG(outputBuffer* buffer, const char * name, ITG * data, unsigned length);
unsigned BufferGetLengthDouble(outputBuffer* buffer, const char * name);
unsigned BufferGetLengthITG(outputBuffer* buffer, const char * name);
void BufferLoadDouble(outputBuffer* buffer, const char * name, double * data, unsigned length);
void BufferLoadITG(outputBuffer* buffer, const char * name, ITG * data, unsigned length);
void BufferNextIter(outputBuffer* buffer);
void BufferClear(outputBuffer* buffer);

#endif
