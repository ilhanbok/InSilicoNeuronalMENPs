% parseRepr2.m
% Description: Parse the data within the files from the LFPy output
% Usage: Place in the same folder as CategoricalPlotGenerator.m
% Author(s): Ilhan Bok
% Last Modified: Feb. 26, 2022

% Parse and store variables from a given magnetization quantification
% produced from Python function `repr`

% Index choices
% 0 = Original Parameters
% 1 = Thick Tufted Variant (L5 TTPC2)
% 2 = Pyramidal Variant (L23 PC)
choices = [  
1
];

% Select input file direction
dir = 'X'

% Ilhan Bok, 2021
prefix1 = 'PYR_VAR_QuantifyMagDist(';
prefix2 = 'PYR_VAR_QuantifyMagDist(';
names = {'H17.06.005.12.12.03_603514493_m'
'H17.06.005.12.15.05_675855243_m'
'H17.06.006.11.10.01_712559850_m'
'H16.03.010.13.05.03_712796148_m'
'H17.06.009.11.04.07_680043920_m'
'H17.03.009.11.11.09_667294513_m'
'H17.06.003.11.03.02_700233253_m'
'H16.06.009.01.01.04.02_595952514_m'
'H16.06.007.01.07.02_572837189_m'
'H16.06.007.01.05.03_571972079_m'
'H16.06.007.01.07.03_709264598_m'
'H16.06.007.01.07.04_582653938_m'
'H17.06.007.11.05.04_603515592_m'
'H16.06.012.11.08.03_666145381_m'
'H16.06.012.11.08.05_572836836_m'
'H16.06.012.11.08.04_668616935_m'
};
suffix = strcat(')',dir,'.out');

for k = 1 : length(names)
  fid = fopen(strcat(prefix1,names{k},suffix));
  raw = fread(fid,inf);
  str=char(raw');
  fclose(fid);
  newStr = strrep(str,'''','"');
  eval(['cell' num2str(k) 'b =  jsondecode(newStr);']);
  fopen(strcat(prefix2,names{k},suffix));
  raw = fread(fid,inf);
  str=char(raw');
  fclose(fid);
  newStr = strrep(str,'''','"');
  eval(['cell' num2str(k) ' =  jsondecode(newStr);']);
end