# Convert line endings from OSX Excel to Linux.
# This command replaces the original file (iin place editing, '-i', option to perl)
# '-P' option iterates over all the input files.
# DJP-S 12-March-2015

perl -pi -e 's/\r//g' $1 
