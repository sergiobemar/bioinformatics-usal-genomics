#!/bin/bash

##############################################################################
# Return the current date tiem in string format without spaces, for namefiles
##############################################################################
function_get_now() {
    date +"%Y-%m-%d %H:%M:%S"
}

##############################################################################
# Return the current date tiem in string format without spaces, for namefiles
##############################################################################
function_get_now_logs() {
    date +"%Y%m%d_%H%M%S"
}

##############################################################################
# Rename all the files of the input directory using lowercase
#
# Arguments:
#   - Input directory
##############################################################################
rename_to_lower_files() {

    echo ""
    echo "$(function_get_now) rename_to_lower_files(): Renaming to lower the files in directory ${1}"
    for i in $( ls ${1} | grep [A-Z] )
    do
        mv -i $i `echo $i | tr 'A-Z' 'a-z'`
    done
}