@echo off

:: loopMorphoL23.bat
:: Description: Runs a biophysical simulation on Allen Brain Atlas morphologies to
::              compute surrounding voxel-based magnetization.
:: Usage: Replace filenames and run from the Windows command line
:: Author(s): Ilhan Bok
:: Last Modified: Jan. 19, 2022

:: MTG3
py .\QuantifyMag_MPIvoxel_L23.py zzzzH17.06.005.12.12.03_603514493_m
py .\QuantifyMag_MPIvoxel_L23.py H17.06.005.12.15.05_675855243_m
py .\QuantifyMag_MPIvoxel_L23.py H17.06.006.11.10.01_712559850_m
py .\QuantifyMag_MPIvoxel_L23.py H16.03.010.13.05.03_712796148_m
:: MTG6
py .\QuantifyMag_MPIvoxel_L23.py H17.06.009.11.04.07_680043920_m
py .\QuantifyMag_MPIvoxel_L23.py H17.03.009.11.11.09_667294513_m
py .\QuantifyMag_MPIvoxel_L23.py H17.06.003.11.03.02_700233253_m
py .\QuantifyMag_MPIvoxel_L23.py H16.06.009.01.01.04.02_595952514_m
:: FL3
py .\QuantifyMag_MPIvoxel_L23.py H16.06.007.01.07.02_572837189_m
py .\QuantifyMag_MPIvoxel_L23.py H16.06.007.01.05.03_571972079_m
py .\QuantifyMag_MPIvoxel_L23.py H16.06.007.01.07.03_709264598_m
py .\QuantifyMag_MPIvoxel_L23.py H16.06.007.01.07.04_582653938_m
:: MFG3
py .\QuantifyMag_MPIvoxel_L23.py H17.06.007.11.05.04_603515592_m
py .\QuantifyMag_MPIvoxel_L23.py H16.06.012.11.08.03_666145381_m
py .\QuantifyMag_MPIvoxel_L23.py H16.06.012.11.08.05_572836836_m
py .\QuantifyMag_MPIvoxel_L23.py H16.06.012.11.08.04_668616935_m