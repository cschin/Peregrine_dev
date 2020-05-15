#!/bin/bash
git checkout 16df8cfe868d5c1598df32d31e250c206618c874
pushd py
rm -rf .eggs/ dist/ build/ peregrine.egg-info/ peregrine_pypy.egg-info get-pip.py
popd

pushd src
make clean
popd

tar czvf src.tgz src/ ksw2/ falcon/ py/ .git/
mv src.tgz docker_paper/

pushd docker_paper/
docker build . --tag cschin/peregrine:hg002_mhc
popd
