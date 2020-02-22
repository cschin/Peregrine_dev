#!/bin/bash

function install_deps() {
    echo "Install deps ..."
    if [[ "$OSTYPE" == "linux-gnu" ]]; then
        # Linux
        sudo apt-get install -y parallel tree
    elif [[ "$OSTYPE" == "darwin"* ]]; then
        # Mac OSX
        brew install parallel tree
    else
        echo "Cannot install deps for $OSTYPE"
    fi
}

if ! hash parallel tree; then
    install_deps
fi

. $CONDAROOT/bin/activate
conda create -n peregrine -y python=3.7

conda activate peregrine

pushd py
rm -rf .eggs/ dist/ build/ peregrine.egg-info/ peregrine_pypy.egg-info get-pip.py
python3 setup.py install
python3 setup.py clean --all
popd
git clone -b peregrine https://github.com/cschin/pypeFLOW.git
pushd pypeFLOW
python3 setup.py install
popd

conda install -c conda-forge -y pypy3.6
conda install -c bioconda -y pysam ncls
pushd py
wget -q https://bootstrap.pypa.io/get-pip.py
pypy3 get-pip.py
pypy3 setup_pypy.py install
popd

pushd src
make all
make install
popd

#python3 -m pip install cffi==1.12.2
