export ALGENCAN=~/algencan-3.1.1

rm -f seropositives

gfortran -O3 -w -fcheck=all -g seropositives.f90 -L$ALGENCAN/lib -lalgencan -lhsl sort.o subset.o -o seropositives

./seropositives
