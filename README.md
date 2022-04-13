# ATC_3D
Source code of DGMM Paper: Tangential cover for 3D irregular noisy digital curves

## Compilation

cd Sources

mkdir build

cd build

cmake .. -DDGtal_DIR=DGTAL_DIR [-DCMAKE_BUILD_TYPE=Release / Debug]

make -j4

## Execution
./ATC3D -i ../Samples/ball.pgm --a11 1.5 --a12 0.2 --a21 0.5 --a22 1.2


## Help
./ATC3D -h

 Adaptive Tangential Cover for 3D curve

 Example:
 
