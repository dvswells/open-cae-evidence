#!/usr/bin/env bash
# Official headless Elmer 9.0, GPL/LGPL; no downloaded executable binaries.
set -euo pipefail
prefix=$(realpath -m "$1")
work=$(realpath -m "$2")
commit=9352634829dfde24439000571d64438d8c683264
mkdir -p "$work" "$prefix"
if ! test -d "$work/source/.git"; then
  git init "$work/source"
  git -C "$work/source" remote add origin https://github.com/ElmerCSC/elmerfem.git
fi
git -C "$work/source" fetch --depth 1 origin "$commit"
git -C "$work/source" checkout --detach FETCH_HEAD
test "$(git -C "$work/source" rev-parse HEAD)" = "$commit"
git -C "$work/source" diff --exit-code
git -C "$work/source" show -s --format=fuller HEAD > "$work/source-provenance.txt"
gcc --version > "$work/compiler.txt"
gfortran --version >> "$work/compiler.txt"
cmake --version >> "$work/compiler.txt"
cmake -S "$work/source" -B "$work/build" \
  -DCMAKE_INSTALL_PREFIX="$prefix" -DCMAKE_BUILD_TYPE=Release \
  -DWITH_MPI=OFF -DWITH_OpenMP=OFF -DWITH_ELMERGUI=OFF \
  -DWITH_ELMERPOST=OFF -DWITH_ElmerIce=OFF -DWITH_CONTRIB=OFF \
  -DWITH_Mumps=OFF -DWITH_Hypre=OFF -DWITH_LUA=OFF \
  2>&1 | tee "$work/configure.log"
cmake --build "$work/build" --parallel 2 --target Solver_TGT ElmerGrid HeatSolve FluxSolver ResultOutputSolve 2>&1 | tee "$work/build.log"
# Install only the required unmodified upstream targets. Upstream's global
# install includes every optional physics module and cannot install a subset.
mkdir -p "$prefix/bin" "$prefix/lib" "$prefix/share/elmersolver/lib"
cp "$work/build/fem/src/ElmerSolver" "$work/build/elmergrid/src/ElmerGrid" "$prefix/bin/"
cp "$work/build/fem/src/"lib{elmersolver,mpi_stubs}.so "$prefix/lib/"
cp "$work/build/matc/src/libmatc.so" "$work/build/fhutiter/src/libfhuti.so" "$work/build/mathlibs/src/arpack/libarpack.so" "$prefix/lib/"
cp "$work/build/fem/src/modules/"{HeatSolve,FluxSolver,ResultOutputSolve}.so "$prefix/share/elmersolver/lib/"
cp "$work/source/fem/src/"{elements.def,SOLVER.KEYWORDS} "$prefix/share/elmersolver/lib/"
# Remove build-tree RPATHs; runtime uses explicit ELMER_HOME and LD_LIBRARY_PATH.
find "$prefix/bin" "$prefix/lib" "$prefix/share/elmersolver/lib" -name '*.so' -o -name ElmerSolver -o -name ElmerGrid | while read -r binary; do
  printf 'file(RPATH_REMOVE FILE "%s")\n' "$binary" >> "$work/remove-rpath.cmake"
done
cmake -P "$work/remove-rpath.cmake" 2>&1 | tee "$work/install.log"
export LD_LIBRARY_PATH="$prefix/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
export ELMER_HOME="$prefix"
cp "$work/build/CMakeCache.txt" "$work/cmake-cache.txt"
cp -r "$work/source/license_texts" "$prefix/source-licenses"
printf '%s\n' "$commit" > "$prefix/source-commit.txt"
find "$prefix" -type f -print0 | sort -z | xargs -0 sha256sum > "$work/binary-manifest.txt"
ldd "$prefix/bin/ElmerSolver" > "$work/solver-libraries.txt"
ldd "$prefix/bin/ElmerGrid" > "$work/grid-libraries.txt"
