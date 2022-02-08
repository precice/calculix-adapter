#!/usr/bin/bash

DISTRIBUTION="_$1"
echo $DISTRIBUTION
ADAPTER_VERSION="2.19.0"
PACKAGE_VERSION="1"

# Compress the changelog, strip the binaries
cp changelog.Debian calculix-precice2_$(ADAPTER_VERSION)-$(PACKAGE_VERSION)_amd64/usr/share/doc/calculix-precice2/changelog.Debian
# Options : --best for best compression, -f for removing file if it was there, -n for no time stamp
gzip --best -f -n calculix-precice2_$(ADAPTER_VERSION)-$(PACKAGE_VERSION)_amd64/usr/share/doc/calculix-precice2/changelog.Debian
strip --strip-unneeded calculix-precice2_$(ADAPTER_VERSION)-$(PACKAGE_VERSION)_amd64/usr/bin/ccx_preCICE

#Compile and compress the manual

pandoc manpage.md -s -t man -o ccx_preCICE.1
mkdir -p calculix-precice2_$(ADAPTER_VERSION)-$(PACKAGE_VERSION)_amd64/usr/share/man/man1
gzip -9 -n -f ccx_preCICE.1
chmod 644 ccx_preCICE.1.gz
mv ccx_preCICE.1.gz calculix-precice2_$(ADAPTER_VERSION)-$(PACKAGE_VERSION)_amd64/usr/share/man/man1

# Copy to a folder with appropriate postfix

cp -r calculix-precice2_$(ADAPTER_VERSION)-$(PACKAGE_VERSION)_amd64/ "calculix-precice2_$(ADAPTER_VERSION)-$(PACKAGE_VERSION)_amd64$DISTRIBUTION"

dpkg-deb --build --root-owner-group "calculix-precice2_$(ADAPTER_VERSION)-$(PACKAGE_VERSION)_amd64$DISTRIBUTION"
lintian *.deb

