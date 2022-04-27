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
   
   ``./ATC3D -i ../data/astroid_noise.dat -m ../MeaningfulThickness/``

<table cellpadding="3">
		<tr>
		<td align="center" valign="center">
			<a href="https://github.com/ngophuc/ATC_3D/blob/main/data/Sinus_ATC3D.png">
				<img height="250" width="450" src="https://github.com/ngophuc/ATC_3D/blob/main/data/Sinus_ATC3D.png" alt="sinus_noise" />
			</a>	
		<br />
		sinus_noise.dat
		</td>
		<td align="center" valign="center">
			<a href="https://github.com/ngophuc/ATC_3D/blob/main/data/Vasque_ATC3D.png">
				<img height="250" width="300" src="https://github.com/ngophuc/ATC_3D/blob/main/data/Vasque_ATC3D.png" alt="vasque_noise" />
			</a>
		<br />
		vasque_noise.dat
		</td>	
		<td align="center" valign="center">
			<a href="https://github.com/ngophuc/ATC_3D/blob/main/data/Astroid_ATC3D.png">
				<img height="250" width="300" src="https://github.com/ngophuc/ATC_3D/blob/main/data/Astroid_ATC3D.png" alt="Astroid_ATC3D" />
			</a>
		<br />
		astroid_noise.dat
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
 
