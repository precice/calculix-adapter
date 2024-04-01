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

/**
 * @brief Get a handle to a buffer object, valid in C parts of the code
 *
 * @return pointer to a buffer
 */
outputBuffer *BufferCreate();

/**
 * @brief Free the buffer object
 *
 * @param buffer
 */
void BufferFree(outputBuffer *buffer);

/**
 * @brief Saves an array of "double" at the current step, under a given name
 *
 * @param buffer handle to the buffer
 * @param name name of the data to save
 * @param data array of values
 * @param length size of the array
 */
void BufferSaveDouble(outputBuffer *buffer, const char *name, double *data, unsigned length);

/**
 * @brief Saves an array of "ITG" at the current step, under a given name
 *
 * @param buffer handle to the buffer
 * @param name name of the data to save
 * @param data array of values
 * @param length size of the array
 */
void BufferSaveITG(outputBuffer *buffer, const char *name, ITG *data, unsigned length);

/**
 * @brief Loads an array of "double" at the current step (modified through ReadNext and Clear) with a given name
 *
 * @param buffer handle to the buffer
 * @param name name of the data to save
 * @param data array of values
 *
 * @pre data is long enough
 */
void BufferLoadDouble(outputBuffer *buffer, const char *name, double *data);

/**
 * @brief Loads an array of "ITG" at the current step (modified through ReadNext and Clear) with a given name
 *
 * @param buffer handle to the buffer
 * @param name name of the data to save
 * @param data array of values
 *
 * @pre data is long enough
 */
void BufferLoadITG(outputBuffer *buffer, const char *name, ITG *data);

/**
 * @brief Returns true (1) if BufferLoadX can be called (i.e. if there is still a step to read)
 *
 */
int BufferCanRead(outputBuffer *buffer);
/**
 * @brief Increment the step counter of data to read
 *
 * @param buffer
 */
void BufferReadNext(outputBuffer *buffer);
/**
 * @brief Start to store a new step
 *
 * @param buffer
 */
void BufferWriteNewStep(outputBuffer *buffer);

/**
 * @brief Clear the buffer and reset the reading counter
 *
 * @param buffer
 */
void BufferClear(outputBuffer *buffer);

/**
 * @brief Check how many states are stored
 *
 * @param buffer
 * @return unsigned number of states that can be read
 */
unsigned BufferStoredStates(outputBuffer *buffer);

#endif
