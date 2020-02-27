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

void ConfigReader_Read(char const * configFilename, char const * participantName, char ** preciceConfigFilename, InterfaceConfig ** interfaces, int * numInterfaces);


#endif
