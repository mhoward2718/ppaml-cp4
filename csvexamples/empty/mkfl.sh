#!/bin/bash -eu

fname=problem-$1-query-$2-metric-$3.csv

if [ $# = 5 ]; then
mv $5 $fname;
fi;

touch $fname
tmpname=$(mktemp -p .)

echo $4 | cat - $fname > $tmpname && mv $tmpname $fname

if [ $# = 5 ]; then
emacs -nw -Q $fname;
fi;
