export ALGENCAN=~/algencan-3.1.1

rm -f seropositives_outliers

gfortran -O3 -w -fcheck=all -g seropositives_outliers.f90 -L$ALGENCAN/lib -lalgencan -lhsl sort.o subset.o -o seropositives_outliers

./seropositives_outliers
