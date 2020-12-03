/**********************************************************************************************
 *                                                                                            *
 *       CalculiX adapter for heat transfer coupling and mechanical FSI using preCICE         *
 *       Heat transfer adapter developed by Lucía Cheung with the support of SimScale GmbH    *
 *                                                                                            *
 *       Adapter extended to fluid-structure interaction by Alexander Rusch                   *
 *                                                                                            *
 *********************************************************************************************/

#ifndef CONFIGREADER_H
#define CONFIGREADER_H

typedef struct InterfaceConfig {
	char * facesMeshName;
	char * nodesMeshName;
	char * patchName;
	int map;
	int numWriteData;
	int numReadData;
	char ** writeDataNames;
	char ** readDataNames;
} InterfaceConfig;

typedef struct AdapterConfig {
  int numInterfaces;
  InterfaceConfig * interfaces;
  char * preciceConfigFilename;
} AdapterConfig;

/**  Reads the Adapter Config
 *
 * @precondition the adapterConfig is a pointer to an uninitialized struct of the type adapterConfig
 *
 * @param[in] configFilename the filename of the adapter config
 * @param[in] participantName the name of the participant
 * @param[inout] adapterConfig a pointer to write the configuration to
 *
 */
void ConfigReader_Read(char const * configFilename, char const * participantName, AdapterConfig * adapterConfig);

/** Frees all internal data held by an adapterConfig
 *
 * @precondition adapterConfig points to an initialized instance of adapterConfig.
 * @precondition ConfigReader_Read was called on adapterConfig.
 * @postcondition all memory held by the struct adapterConfig is freed
 * @note This function does not free the pointer adapterConfig
 */
void AdapterConfig_Free(AdapterConfig * adapterConfig);

#endif
