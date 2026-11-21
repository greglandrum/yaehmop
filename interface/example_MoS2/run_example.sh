#!/bin/bash
# Full example: YAeHMOP (bind) -> yaehmop2wtb.py -> WanTiBEXOS (wtb.x) -> comparison.
# Run it inside interface/example_MoS2.  The programs are taken from the PATH, or from
#   BIND=/path/to/bind  WTB=/path/to/wtb.x  ./run_example.sh
# ("bind" is also a bash builtin, so the binary is looked up with type -P.)
set -e
BIND=${BIND:-$(type -P bind)}
WTB=${WTB:-wtb.x}
[ -x "$BIND" ] || { echo "bind not found: set BIND=/path/to/bind"; exit 1; }
ulimit -s unlimited || true
export OMP_NUM_THREADS=1 KMP_STACKSIZE=500mb
here=$(pwd)

# 1) YAeHMOP runs (with and without SOC): band path + H(k), S(k) on the 11x11x1 mesh
for case in soc nosoc; do
  mkdir -p yaehmop_$case && cp MoS2.bind yaehmop_$case/
  cd yaehmop_$case
  if [ $case = soc ]; then
    "$BIND" --soc ../MoS2_soc.dat MoS2.bind ../MoS2_params.dat > /dev/null
  else
    "$BIND" MoS2.bind ../MoS2_params.dat > /dev/null
  fi
  cd "$here"
done

# 2) conversion to the WanTiBEXOS format and 3) wtb.x band structure
for case in soc nosoc; do
  python3 ../yaehmop2wtb.py yaehmop_$case/MoS2.bind -o wtb_$case
  cd wtb_$case
  mkdir -p out
  kind=$(ls input-*.dat | sed 's/input-//; s/.dat//')
  "$WTB" < input-$kind.dat > wtb.log
  cd "$here"
done

# 4) comparison of the bands and figure
python3 compare_bands.py
