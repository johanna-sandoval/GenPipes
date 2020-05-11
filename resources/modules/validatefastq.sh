#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

SOFTWARE=validatefastq
VERSION=0.1.1
JAR=${SOFTWARE}-assembly-${VERSION}.jar
ARCHIVE=$JAR
SOFTWARE_DIR=${SOFTWARE}-${VERSION}
ARCHIVE_URL=https://github.com/biopet/${SOFTWARE}/releases/download/v${VERSION}/$JAR

build() {
  cd $INSTALL_DOWNLOAD
  mkdir -p $INSTALL_DIR/$SOFTWARE_DIR
  cp $JAR $INSTALL_DIR/$SOFTWARE_DIR/ 
}

module_file() {
echo "\
#%Module1.0
proc ModulesHelp { } {
  puts stderr \"\tMUGQIC - $SOFTWARE \"
}
module-whatis \"$SOFTWARE\"

set             root                $INSTALL_DIR/$SOFTWARE_DIR
prepend-path    PATH	            \$root
setenv          VALIDATEFASTQ_JAR   \$root/$JAR
"
}

# Call generic module install script once all variables and functions have been set
MODULE_INSTALL_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $MODULE_INSTALL_SCRIPT_DIR/install_module.sh $@