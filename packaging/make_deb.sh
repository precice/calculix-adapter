#!/usr/bin/bash
set -e -u

DISTRIBUTION="_$1"
ADAPTER_VERSION="2.20.2"
PACKAGE_VERSION="1"

PACKAGE_FOLDER="calculix-precice3_$ADAPTER_VERSION-${PACKAGE_VERSION}_amd64"

# Compress the changelog, strip the binaries
cp changelog.Debian debian/usr/share/doc/calculix-precice3/changelog.Debian
# Options : --best for best compression, -f for removing file if it was there, -n for no time stamp
gzip --best -f -n debian/usr/share/doc/calculix-precice3/changelog.Debian
strip --strip-unneeded debian/usr/bin/ccx_preCICE

#Compile and compress the manual

pandoc manpage.md -s -t man -o ccx_preCICE.1
mkdir -p debian/usr/share/man/man1
gzip -9 -n -f ccx_preCICE.1
chmod 644 ccx_preCICE.1.gz
mv ccx_preCICE.1.gz debian/usr/share/man/man1

# Copy to a folder with appropriate postfix

cp -r debian/ "$PACKAGE_FOLDER$DISTRIBUTION"

dpkg-deb --build --root-owner-group "$PACKAGE_FOLDER$DISTRIBUTION"
lintian ./*.deb

