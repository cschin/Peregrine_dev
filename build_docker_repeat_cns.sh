#!/bin/bash
pushd py
rm -rf .eggs/ dist/ build/ peregrine.egg-info/ peregrine_pypy.egg-info get-pip.py
popd

pushd src
make clean
popd

#tar czvf src.tgz src/ ksw2/ falcon/ py/ .git/
#mv src.tgz docker_repeat_cns/

pushd docker_repeat_cns/
docker build . --tag cschin/peregrine:repeat_cns
popd
