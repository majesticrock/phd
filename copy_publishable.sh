#!/usr/bin/env bash
set -euo pipefail

mkdir -p ../publishable
mkdir -p ../publishable/cpp

# PhdUtility
SRC=./PhdUtility
DEST=../../publishable/PhdUtility

cd "$SRC"
git ls-files \
    -- . \
    ':(exclude).*' \
    ':(exclude)**/docs/**' \
| rsync -a --files-from=- "./" "$DEST/"


########
cd ../cpp

# FermionCommute
SRC=./FermionCommute
DEST=../../../publishable/cpp/FermionCommute

cd "$SRC"
git ls-files \
    -- . \
    ':(exclude).*' \
| rsync -a --files-from=- "./" "$DEST/"
cd ..

# ContinuumSystem
SRC=./ContinuumSystem
DEST=../../../publishable/cpp/ContinuumSystem

cd "$SRC"
git ls-files \
    -- . \
    ':(exclude).*' \
    ':(exclude)*.sh' \
    ':(exclude)*.slurm' \
    ':(exclude)**/benchmark.cpp' \
    ':(exclude)**/NoScreeningIntegrals.cpp' \
| rsync -a --files-from=- "./" "$DEST/"
cd ..

# Hubbard
SRC=./Hubbard
DEST=../../../publishable/cpp/Hubbard

cd "$SRC"
git ls-files \
    -- . \
    ':(exclude).*' \
    ':(exclude)*.sh' \
    ':(exclude)*.slurm' \
| rsync -a --files-from=- "./" "$DEST/"
cd ..

# LatticeCUT
SRC=./LatticeCUT
DEST=../../../publishable/cpp/LatticeCUT

cd "$SRC"
git ls-files \
    -- . \
    ':(exclude).*' \
    ':(exclude)*.sh' \
    ':(exclude)*.slurm' \
| rsync -a --files-from=- "./" "$DEST/"
cd ../..

########

# superproject

DEST=../publishable

rsync -a \
    --include='plot_examples/***' \
    --include='README.md' \
    --include='TEST_RESULTS.md' \
    --include='test_everything.sh' \
    --exclude='*' \
    "./" "$DEST/"