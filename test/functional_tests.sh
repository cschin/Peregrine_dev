#!/bin/bash

test -e ssshtest || wget -q https://raw.githubusercontent.com/ryanlayer/ssshtest/master/ssshtest

. ssshtest
set -uo pipefail

# ecoli_K12
pushd ./test/ecoli_K12

echo "[TEST] ecoli_K12 simreads ..."
run ecoli_K12_simreads bash -c "make simreads && echo 'SUCCESS'"
assert_in_stdout "SUCCESS"

echo "[TEST] ecoli_K12 test ..."
run ecoli_K12_test bash -c "make test && echo 'SUCCESS'"
assert_in_stdout "SUCCESS"

# Just get a visual validation of output
cat logs/all.log
ls -tal wd/asm

popd
