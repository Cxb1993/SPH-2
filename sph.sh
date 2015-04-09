#! /bin/bash

pathname="$(dirname "$(readlink -f "$0")")"

for ((n=0;n<$1;n++))
do
    $pathname/sph
done
