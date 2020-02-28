/**********************************************************************************************
 *                                                                                            *
 *       CalculiX adapter for heat transfer coupling and mechanical FSI using preCICE         *
 *       Heat transfer adapter developed by Lucía Cheung with the support of SimScale GmbH    *
 *                                                                                            *
 *       Adapter extended to fluid-structure interaction by Alexander Rusch                   *
 *                                                                                            *
 *********************************************************************************************/

#include "ConfigReader.hpp"
#include "yaml-cpp/yaml.h"
#include <iostream>
#include <cstring>
#include <algorithm>
#include <cstdlib>
#include <iterator>

void ConfigReader_Read( char const * configFilename, char const * participantName, char ** preciceConfigFilename, InterfaceConfig ** interfaces, int * numInterface )
{

	YAML::Node config = YAML::LoadFile( configFilename );

	*preciceConfigFilename = strdup( config["precice-config-file"].as<std::string>().c_str() );

	*numInterface = config["participants"][participantName]["interfaces"].size();
	*interfaces = (InterfaceConfig*) malloc( *numInterface * sizeof( InterfaceConfig ) );

	

	for( int i = 0 ; i < *numInterface ; i++ )
	{
    InterfaceConfig * currentInterfacePointer = *interfaces;
    std::advance( currentInterfacePointer, i );
    new ( currentInterfacePointer ) InterfaceConfig();
    InterfaceConfig& interface = *currentInterfacePointer;

		if( config["participants"][participantName]["interfaces"][i]["nodes-mesh"] )
		{
			interface.nodesMeshName = strdup( config["participants"][participantName]["interfaces"][i]["nodes-mesh"].as<std::string>().c_str() );
		}
		else if ( config["participants"][participantName]["interfaces"][i]["nodes-mesh-with-connectivity"] ) 
		{
			interface.nodesMeshName = strdup( config["participants"][participantName]["interfaces"][i]["nodes-mesh-with-connectivity"].as<std::string>().c_str() );
			interface.map = 1;
		}
		else
		{
			interface.nodesMeshName = NULL;
		}

		if( config["participants"][participantName]["interfaces"][i]["faces-mesh"] )
		{
			interface.facesMeshName = strdup( config["participants"][participantName]["interfaces"][i]["faces-mesh"].as<std::string>().c_str() );
		}
		else
		{
			interface.facesMeshName = NULL;
		}
		
		if( config["participants"][participantName]["interfaces"][i]["mesh"] )
		{
			interface.facesMeshName = strdup( config["participants"][participantName]["interfaces"][i]["mesh"].as<std::string>().c_str() );
		}

		std::string patchName = config["participants"][participantName]["interfaces"][i]["patch"].as<std::string>();
		std::transform( patchName.begin(), patchName.end(), patchName.begin(), toupper );
		interface.patchName = strdup( patchName.c_str() );
		
		interface.numWriteData = config["participants"][participantName]["interfaces"][i]["write-data"].size();
		interface.numReadData = config["participants"][participantName]["interfaces"][i]["read-data"].size();

		if( interface.numWriteData == 0 )
		{
			// write-data is a string
			interface.numWriteData = 1;
			interface.writeDataNames = (char**) malloc( sizeof( char* ) * interface.numWriteData );
			interface.writeDataNames[0] = strdup( config["participants"][participantName]["interfaces"][i]["write-data"].as<std::string>().c_str() );
		}
		else
		{
			// write-data is an array
			interface.writeDataNames = (char**) malloc( sizeof( char* ) * interface.numWriteData );

			for( int j = 0 ; j < interface.numWriteData ; j++ )
			{
				interface.writeDataNames[j] = strdup( config["participants"][participantName]["interfaces"][i]["write-data"][j].as<std::string>().c_str() );
			}
		}

		if( interface.numReadData == 0 )
		{
			// read-data is a string
			interface.numReadData = 1;
			interface.readDataNames = (char**) malloc( sizeof( char* ) * interface.numReadData );
			interface.readDataNames[0] = strdup( config["participants"][participantName]["interfaces"][i]["read-data"].as<std::string>().c_str() );
		}
		else
		{
			// read-data is an array
			interface.readDataNames = (char**) malloc( sizeof( char* ) * interface.numReadData );

			for( int j = 0 ; j < interface.numReadData ; j++ )
			{
				interface.readDataNames[j] = strdup( config["participants"][participantName]["interfaces"][i]["read-data"][j].as<std::string>().c_str() );
			}
		}	

	}
}

void InterfaceConfig_Free(InterfaceConfig * interface)
{
  if (interface == nullptr) return;

  // Owning char arrays
	free(interface->facesMeshName);
	free(interface->nodesMeshName);
	free(interface->patchName);

  // Owning arrays of char arrays
  char ** writeDataNamesEnd = interface->writeDataNames;
  std::advance(writeDataNamesEnd, interface->numWriteData);
  std::for_each(interface->writeDataNames, writeDataNamesEnd, free);
  free(interface->writeDataNames);

  char ** readDataNamesEnd = interface->readDataNames;
  std::advance(readDataNamesEnd, interface->numReadData);
  std::for_each(interface->readDataNames, readDataNamesEnd, free);
  free(interface->readDataNames);

  // The interface itself
  free(interface);
}
