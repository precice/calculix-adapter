#!/usr/bin/bash

DISTRIBUTION="_$1"
ADAPTER_VERSION="2.19.0"
PACKAGE_VERSION="1"

# First folder is the one where we put data, second is created when building the package
DEFAULT_FOLDER="calculix-precice2_pkg"
PACKAGE_FOLDER="calculix-precice2_$ADAPTER_VERSION-${PACKAGE_VERSION}_amd64"

# Compress the changelog, strip the binaries
cp changelog.Debian $DEFAULT_FOLDER/usr/share/doc/calculix-precice2/changelog.Debian
# Options : --best for best compression, -f for removing file if it was there, -n for no time stamp
gzip --best -f -n $DEFAULT_FOLDER/usr/share/doc/calculix-precice2/changelog.Debian
strip --strip-unneeded $DEFAULT_FOLDER/usr/bin/ccx_preCICE

#Compile and compress the manual

pandoc manpage.md -s -t man -o ccx_preCICE.1
mkdir -p $DEFAULT_FOLDER/usr/share/man/man1
gzip -9 -n -f ccx_preCICE.1
chmod 644 ccx_preCICE.1.gz
mv ccx_preCICE.1.gz $DEFAULT_FOLDER/usr/share/man/man1

# Copy to a folder with appropriate postfix

cp -r $DEFAULT_FOLDER/ "$PACKAGE_FOLDER$DISTRIBUTION"

dpkg-deb --build --root-owner-group "$PACKAGE_FOLDER$DISTRIBUTION"
lintian ./*.deb

