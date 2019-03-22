The frdToVTKConverter converst a Calculix results output file (in .frd format) to a VTK file that can be read by Paraview. The converter was originally published on FreeCAD forum by user mamah:

https://forum.freecadweb.org/viewtopic.php?t=11664

The original converter was written in C# split over various files, and has been adapted into a single executable. The software Mono must be downloaded, which allows the running of C# executables in Ubuntu. The following commands must be run in a terminal to ensure the conversion of FRD to VTK. Disclosure: This converter has only been tested on Ubuntu 18.04.

sudo apt update
sudo apt install mono-complete

To make an executable of the supplied file "frdToVTKConverter.cs", run:

mcs -out:frdToVTKConverter.exe frdToVTKConverter.cs

Then to run the converter to produce VTK files from an FRD file, run:

mono frdToVTKConverter.exe /path/to/frdFile.frd

When performing a conversion, the adapter will search for all data types that was saved in the .frd file (e.g. displacement, stress, strain etc.). After reading this list it will pause with the total number of timesteps that were read. The user must re-enter the number of timesteps into the terminal for the writing of the VTK files.
