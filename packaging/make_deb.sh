#!/usr/bin/bash

DISTRIBUTION="_$1"
echo $DISTRIBUTION

# Compress the changelog, strip the binaries
cp changelog.Debian calculix-precice_2.17-1_amd64/usr/share/doc/calculix-precice/changelog.Debian
# Options : -f for removing file if it was there, -n for no time stamp
gzip -9 calculix-precice_2.17-1_amd64/usr/share/doc/calculix-precice/changelog.Debian -f -n
strip --strip-unneeded calculix-precice_2.17-1_amd64/usr/bin/ccx_preCICE

#Compile and compress the manual

pandoc manpage.md -s -t man -o ccx_preCICE.1
mkdir -p calculix-precice_2.17-1_amd64/usr/share/man/man1
gzip -9 -n -f ccx_preCICE.1
chmod 644 ccx_preCICE.1.gz
mv ccx_preCICE.1.gz calculix-precice_2.17-1_amd64/usr/share/man/man1

# Copy to a folder with appropriate postfix

cp -r calculix-precice_2.17-1_amd64/ "calculix-precice_2.17-1_amd64$DISTRIBUTION"

dpkg-deb --build --root-owner-group "calculix-precice_2.17-1_amd64$DISTRIBUTION"
lintian *.deb

