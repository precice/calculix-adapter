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

void ConfigReader_Read( char const * configFilename, char const * participantName, char ** preciceConfigFilename, InterfaceConfig ** interfaces, int * numInterface )
{

	YAML::Node config = YAML::LoadFile( configFilename );

	*preciceConfigFilename = strdup( config["precice-config-file"].as<std::string>().c_str() );

	*numInterface = config["participants"][participantName]["interfaces"].size();
	*interfaces = (InterfaceConfig*) malloc( sizeof( InterfaceConfig ) * *numInterface );

	

	for( int i = 0 ; i < *numInterface ; i++ )
	{
		if( config["participants"][participantName]["interfaces"][i]["nodes-mesh"] )
		{
			( *interfaces )[i].nodesMeshName = strdup( config["participants"][participantName]["interfaces"][i]["nodes-mesh"].as<std::string>().c_str() );
		}
		else if ( config["participants"][participantName]["interfaces"][i]["nodes-mesh-with-connectivity"] ) 
		{
			( *interfaces )[i].nodesMeshName = strdup( config["participants"][participantName]["interfaces"][i]["nodes-mesh-with-connectivity"].as<std::string>().c_str() );
			( *interfaces )[i].map = 1;
		}
		else
		{
			( *interfaces )[i].nodesMeshName = NULL;
		}

		if( config["participants"][participantName]["interfaces"][i]["faces-mesh"] )
		{
			( *interfaces )[i].facesMeshName = strdup( config["participants"][participantName]["interfaces"][i]["faces-mesh"].as<std::string>().c_str() );
		}
		else
		{
			( *interfaces )[i].facesMeshName = NULL;
		}
		
		if( config["participants"][participantName]["interfaces"][i]["mesh"] )
		{
			( *interfaces )[i].facesMeshName = strdup( config["participants"][participantName]["interfaces"][i]["mesh"].as<std::string>().c_str() );
		}

		std::string patchName = config["participants"][participantName]["interfaces"][i]["patch"].as<std::string>();
		std::transform( patchName.begin(), patchName.end(), patchName.begin(), toupper );
		( *interfaces )[i].patchName = strdup( patchName.c_str() );
		
		( *interfaces )[i].numWriteData = config["participants"][participantName]["interfaces"][i]["write-data"].size();
		( *interfaces )[i].numReadData = config["participants"][participantName]["interfaces"][i]["read-data"].size();

		if( ( *interfaces )[i].numWriteData == 0 )
		{
			// write-data is a string
			( *interfaces )[i].numWriteData = 1;
			( *interfaces )[i].writeDataNames = (char**) malloc( sizeof( char* ) * ( *interfaces )[i].numWriteData );
			( *interfaces )[i].writeDataNames[0] = strdup( config["participants"][participantName]["interfaces"][i]["write-data"].as<std::string>().c_str() );
		}
		else
		{
			// write-data is an array
			( *interfaces )[i].writeDataNames = (char**) malloc( sizeof( char* ) * ( *interfaces )[i].numWriteData );

			for( int j = 0 ; j < ( *interfaces )[i].numWriteData ; j++ )
			{
				( *interfaces )[i].writeDataNames[j] = strdup( config["participants"][participantName]["interfaces"][i]["write-data"][j].as<std::string>().c_str() );
			}
		}

		if( ( *interfaces )[i].numReadData == 0 )
		{
			// read-data is a string
			( *interfaces )[i].numReadData = 1;
			( *interfaces )[i].readDataNames = (char**) malloc( sizeof( char* ) * ( *interfaces )[i].numReadData );
			( *interfaces )[i].readDataNames[0] = strdup( config["participants"][participantName]["interfaces"][i]["read-data"].as<std::string>().c_str() );
		}
		else
		{
			// read-data is an array
			( *interfaces )[i].readDataNames = (char**) malloc( sizeof( char* ) * ( *interfaces )[i].numReadData );

			for( int j = 0 ; j < ( *interfaces )[i].numReadData ; j++ )
			{
				( *interfaces )[i].readDataNames[j] = strdup( config["participants"][participantName]["interfaces"][i]["read-data"][j].as<std::string>().c_str() );
			}
		}	

	}
}

