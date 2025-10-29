#!/bin/bash
RESOURCES_DIR=${1:-"Resources"};

#---------------------------------------------------------------------------
# Constants
#---------------------------------------------------------------------------
SCRIPT_DIR="$(dirname "$(readlink -f "$0")")"
CMPILEUP_DIR="${RESOURCES_DIR}/Sites"

#---------------------------------------------------------------------------
# Copy file to resources directory
#---------------------------------------------------------------------------
cp ${SCRIPT_DIR}/CMPileup.ini ${CMPILEUP_DIR}/CMPileup.ini
cp ${SCRIPT_DIR}/CMPileup.stranded.ini ${CMPILEUP_DIR}/CMPileup.stranded.ini