#!/usr/bin/env python3
# -*- coding:utf-8 -*-
#
# test_ctg_aln.py
# peregrine
#
# Created by Haibao Tang on 11/27/19
#

import pytest


def get_mmer_index(py_mmer):
    mmer_index = {}
    for i in range(py_mmer.mmers.n):
        mmer = py_mmer.mmers.a[i]
        rid = mmer.y >> 32
        mmer_index.setdefault(rid, [None, None])
        if mmer_index[rid][0] == None:
            mmer_index[rid][0] = i
        if mmer_index[rid][1] == None or i > mmer_index[rid][1]:
            mmer_index[rid][1] = i
    return mmer_index


def test_shimmer4py():
    from peregrine._shimmer4py import ffi as shimmer_ffi

    from peregrine._shimmer4py import lib as shimmer4py

    py_mmer_L2 = shimmer_ffi.new("py_mmer_t *")
    shimmer4py.build_shimmer_map4py(
        py_mmer_L2,
        b"test/ecoli_K12/wd/index/seq_dataset",
        b"test/ecoli_K12/wd/index/shmr-L2",
        1,
        1,
        2,
        240,
    )

    mmer_index_L2 = get_mmer_index(py_mmer_L2)
    assert len(mmer_index_L2) == pytest.approx(5000, 100)
