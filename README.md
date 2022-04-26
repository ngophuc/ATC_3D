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

Tangential cover for 3D irregular noisy digital curves.

Example:
 	 ATC3D --input <FileName> --imaGeneDir <imaGeneDir> 

Usage: ./ATC3D [OPTIONS] 1

Positionals:
  1 TEXT:FILE REQUIRED                  Input file.

Options:
  -h,--help                             Print this help message and exit
  -i,--input TEXT:FILE REQUIRED         Input file.
  -m,--mt TEXT                          MeaningfulThickness directory for noise detection (default ../MeaningfulThickness/).
 
