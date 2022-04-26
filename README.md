# ATC_3D
Source code of DGMM Paper: Tangential cover for 3D irregular noisy digital curves

## Compilation

cd Sources

mkdir build

cd build

cmake .. -DDGtal_DIR=DGTAL_DIR [-DCMAKE_BUILD_TYPE=Release / Debug]

make -j4

## Execution
./ATC3D -i ../data/sinus_noise.dat -m ../MeaningfulThickness/

./ATC3D -i ../data/vasque_noise.dat -m ../MeaningfulThickness/

## Help
./ATC3D -h

 Adaptive Tangential Cover for 3D curve

 Example:
 
