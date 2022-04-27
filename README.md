Source code of DGMM Paper: 
# Tangential cover for 3D irregular noisy digital curves

## Compilation

1. Build meaningful thickness noise detector

 	 ``cd <SourcesDirectory>/MeaningfulThickness; mkdir build; cd build``

 	 ``cmake ..; make -j4``

2. Build 3D adaptive tangential cover 

   ``cd <SourcesDirectory>; mkdir build; cd build``

 	 ``cmake .. -DDGtal_DIR=<DGTAL_DIR> [-DCMAKE_BUILD_TYPE=Release / Debug]``
  
   ``make -j4``

## Execution
   ``./ATC3D -i ../data/sinus_noise.dat -m ../MeaningfulThickness/``

   ``./ATC3D -i ../data/vasque_noise.dat -m ../MeaningfulThickness/``

<table cellpadding="3">
		<tr>
		<td align="center" valign="center">
			<a href="https://github.com/ngophuc/ATC_3D/data/bird5.png">
				<img width="300" src="https://github.com/ngophuc/ModifiedAdaptiveTagentialCover/blob/master/Samples/bird5.png" alt="Input image" />
			</a>	
		<br />
		sinus_noise.dat
		</td>
		<td align="center" valign="center">
			<a href="https://github.com/ngophuc/ModifiedAdaptiveTagentialCover/blob/master/Results/bird5_ATC.pdf">
				<img width="300" src="https://github.com/ngophuc/ModifiedAdaptiveTagentialCover/blob/master/Results/bird5_ATC.png" alt="Adaptive Tagential Cover result" />
			</a>
		<br />
		vasque_noise.dat
		</td>	
		<td align="center" valign="center">
			<a href="https://github.com/ngophuc/ModifiedAdaptiveTagentialCover/blob/master/Results/bird5_DPnew_ATC.pdf">
				<img width="300" src="https://github.com/ngophuc/ModifiedAdaptiveTagentialCover/blob/master/Results/bird5_DPnew_ATC.png" alt="Polygonal approximation with ATC" />
			</a>
		<br />
		vasque_noise.dat
		</td>			
		</tr>
	</table>
 
## Help
``./ATC3D -h``

Tangential cover for 3D irregular noisy digital curves.

Example:

 	 ATC3D --input <FileName> --mt <MeaningfulThicknessDir> 

Usage: ./ATC3D [OPTIONS] 1

Positionals:
  1 TEXT:FILE REQUIRED                  Input file.

Options:

  -h,--help                             Print this help message and exit
  
  -i,--input TEXT:FILE REQUIRED         Input file.
  
  -m,--mt TEXT                          MeaningfulThickness directory for noise detection (default ../MeaningfulThickness/).
 
