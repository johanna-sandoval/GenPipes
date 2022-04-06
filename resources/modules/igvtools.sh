#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

SOFTWARE=igvtools
VERSION=2.11.9
ARCHIVE=${SOFTWARE}_$VERSION.zip
ARCHIVE_URL=http://www.broadinstitute.org/igv/projects/downloads/2.11/$ARCHIVE
SOFTWARE_DIR=$SOFTWARE-$VERSION

build() {
  cd $INSTALL_DOWNLOAD
  unzip $ARCHIVE
  mv IGVTools $INSTALL_DIR/$SOFTWARE_DIR
}

module_file() {
echo "\
#%Module1.0
proc ModulesHelp { } {
  puts stderr \"\tMUGQIC - $SOFTWARE \"
}
module-whatis \"$SOFTWARE\"

set             root                $INSTALL_DIR/$SOFTWARE_DIR
prepend-path    PATH                \$root
setenv          IGVTOOLS_JAR        \$root/igvtools.jar
"
}

# Call generic module install script once all variables and functions have been set
MODULE_INSTALL_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $MODULE_INSTALL_SCRIPT_DIR/install_module.sh $@
