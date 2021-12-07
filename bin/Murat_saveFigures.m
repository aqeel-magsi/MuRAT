function  Murat_saveFigures(figureName,Path)
% function  Murat_plotSections(figure,Path)
% PLOTS sections of the desired model
%
%	Input Parameters:
%       figureName:             Matlab 3D plot
%       Path:                   name of sections (tif) and figure (fig)
%
%   Output:
%       Three sections and 1 3D figure
%
SetFDefaults
saveas(figureName,Path);    
view(90,0)
saveas(figureName,[Path '_SN'], 'tif');
view(0,0)
saveas(figureName,[Path '_WE'], 'tif');
view(0,90)
axis equal
saveas(figureName,[Path '_H'], 'tif');
close(figureName)
end