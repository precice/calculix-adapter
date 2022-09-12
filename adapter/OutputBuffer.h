/*
    Code defining buffering mechanism for the output.

    CalculiX outputs .frd file using the frd() function, which takes many parameters.
    When doing implicit coupling with subcycling in preCICE, only the sub-steps of the final iteration of 
    a time window must be outputted. Thus, it is necessary to store data until the end of the time window:
    - If the window converged, the stored steps can be outputted all at once.
    - If the window must be restarted, the buffer should be cleared.

    All of this is necessary because CalculiX does not expect to go back by more than one step.

    It is assumed that C arrays are large enough: in practice, their size is independent of the step, so this shouldn't
    be an issue.
*/

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

typedef struct outputBuffer outputBuffer;

outputBuffer *BufferCreate();
void          BufferFree(outputBuffer *buffer);
void          BufferSaveDouble(outputBuffer *buffer, const char *name, double *data, unsigned length);
void          BufferSaveITG(outputBuffer *buffer, const char *name, ITG *data, unsigned length);

void BufferLoadDouble(outputBuffer *buffer, const char *name, double *data);
void BufferLoadITG(outputBuffer *buffer, const char *name, ITG *data);

int      BufferCanRead(outputBuffer *buffer);
void     BufferReadNext(outputBuffer *buffer);
void     BufferWriteNewStep(outputBuffer *buffer);
void     BufferClear(outputBuffer *buffer);
unsigned BufferStoredStates(outputBuffer *buffer);

#endif
