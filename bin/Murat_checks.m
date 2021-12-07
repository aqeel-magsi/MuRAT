% ADDITIONAL input variables that are not set by the user.
function Murat                  =   Murat_checks(Murat,mag,timespan)

% INPUTS
dataDirectory                   =   ['./' Murat.input.dataDirectory];
PTime                           =   ['SAChdr.times.' Murat.input.PTime];
PorS                            =   Murat.input.POrS;

origin                          =   Murat.input.origin;
ending                          =   Murat.input.end;
nLat                            =   Murat.input.gridLat;
nLong                           =   Murat.input.gridLong;
nzc                             =   Murat.input.gridZ;
availableVelocity               =   Murat.input.availableVelocity;
velocityModel                   =   ['velocity_models/',Murat.input.namev];

if isempty(Murat.input.originTime)
    originTime                  =   'SAChdr.times.o';
else
    originTime                  =...
        ['SAChdr.times.' Murat.input.originTime];
end

if isempty(Murat.input.STime)
    STime                  =   'SAChdr.times.t0';
else
    STime                  =...
        ['SAChdr.times.' Murat.input.STime];
end

Murat.input.originTime          =   originTime;
Murat.input.PTime               =   PTime;
Murat.input.STime               =   STime;


if exist('./temp','dir')==7    
    delete('./temp/*')
else
    mkdir('./temp')
end

% Checking data
[Murat.input.listSac,~]         =	createsList(Murat,mag,timespan);
[Murat.input.header,flag]       =...
    Murat_testData(dataDirectory,originTime,PTime,STime);

if isequal(flag,1)
    warning('Missing origin times.')
end

if isequal(flag,2)
    warning('Missing S-wave times.')
end
%% VELOCITY MODELS: ORIGINAL, INVERSION, and PROPAGATION
% Save x,y,z in degrees switching as longitude comes second
% Find distance and azimuth to change in meters - requires longitude first
dist_x                          =   deg2km(ending(2)-origin(2))*1000;
dist_y                          =   deg2km(ending(1)-origin(1))*1000;

% Coordinates of inversion points in meters
xM                              =   linspace(0,dist_x,nLong)';
yM                              =   linspace(0,dist_y,nLat)';
zM                              =   linspace(origin(3),ending(3),nzc)';
modvXYZ                         =   Murat_unfoldXYZ(xM,yM,zM);

Murat.input.x                   =   linspace(origin(2),ending(2),nLong)';
Murat.input.y                   =   linspace(origin(1),ending(1),nLat)';
Murat.input.z                   =   linspace(origin(3),ending(3),nzc)';
Murat.input.gridStepX           =   xM(2)-xM(1);
Murat.input.gridStepY           =   yM(2)-yM(1);

modvOriginal                    =   load(velocityModel);
if availableVelocity ==  0

    gridPropagation.x           =   xM';
    gridPropagation.y           =   yM';
    gridPropagation.z           =   zM';
    
    [modv,pvel]                 =...
        Murat_modv1D(modvXYZ,modvOriginal,PorS);
    
    Murat.input.modv            =   modv;
    Murat.input.modvp           =   modv;
    Murat.input.modvPlot        =   [];
    
elseif availableVelocity ==  1
    
    [modvP,modvI,modvIP,pvel]          =...
        Murat_modv3D(modvXYZ,modvOriginal,origin,0);
    
    gridPropagation.x           =   unique(modvP(:,1));
    gridPropagation.y           =   unique(modvP(:,2));
    gridPropagation.z           =   sort(unique(modvP(:,3)),'descend');
    
    Murat.input.modv            =   modvI;
    Murat.input.modvp           =   modvP(:,1:4);    
    Murat.input.modvPlot        =   modvIP;
    
end

Murat.input.gridPropagation     =   gridPropagation;
Murat.input.pvel                =   pvel;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [listWithFolder,listNoFolder]...
    =   createsList(Murat,mag,timespan)
% CREATES a list of visible files in a folder, outputs both with and
% without folder
%
if isempty(mag) && isempty(timespan)
    list                            =   dir(Murat.input.dataDirectory);
    list                            =   list(~startsWith({list.name}, '.'));
    
    listWithFolder                  =	fullfile({list.folder},{list.name})';
    listNoFolder                    =   {list.name}';
    
else
    %% Miriam's option to only read events of certain magnitude
    
    DFolder                     = Murat.input.dataDirectory;
    load('chosen_events.mat')
    orig_data = length(TZ_events.year);
    
    % delete traces outside main area of interest
    to_del = zeros(length(TZ_events.year),1);
    for i=1:length(TZ_events.year)
        if TZ_events.lat(i) > -2.45 || TZ_events.lat(i) < -2.95
            to_del(i) = 1;
        end
        if TZ_events.lon(i) < 35.8  || TZ_events.lon(i) > 36.25
            to_del(i) = 1;
        end
    end
    del_row = find(to_del==1);
    TZ_events(del_row,:)=[];
    disp(['removed ' , num2str(length(find(to_del==1))), ' events because they are far from network'])
    % delete traces below magnitude threshold
    find_mags = find(TZ_events.mag<mag);
    TZ_events(find_mags,:) = [];
    disp(['removed ' , num2str(length(find_mags)), ' events because they are below the magnitude threshold'])
    % keep traces inside chosen time frame
    ex_dates=cellstr(TZ_events.of);
    dates = cellfun(@(x) x(1:6),ex_dates,'un',0);
    find_dates = contains(dates,timespan);
    TZ_events(~find_dates,:) = [];
    disp(['selected ' , num2str(length(find(find_dates==1))), ' events inside the set timeframe'])
    % delete traces below certain depth
    find_depth = find(TZ_events.depth>21);
    TZ_events(find_depth,:) = [];
    disp(['removed ' , num2str(length(find_depth)), ' events because they are below 21 km depth'])
    find_depth = find(TZ_events.depth<1);
    TZ_events(find_depth,:) = [];
    disp(['removed ' , num2str(length(find_depth)), ' events because they are above 1 km depth'])
    % delete events which have only few traces
    find_ns = find(TZ_events.ns<10);
    TZ_events(find_ns,:) = [];
    disp(['removed ' , num2str(length(find_ns)), ' events because they were recorded by less than 10 components'])
    
    disp(['keeping ',num2str(length(TZ_events.lat)), 'of ', num2str(orig_data),' events'])
    
    new_filenames = [];
    new_filedir = [];
    
    for ii = 1:length(TZ_events.year)
        tmp_filename = sprintf('%4d%02d%02d%02d%02d%02d',TZ_events.year(ii),TZ_events.month(ii),...
            TZ_events.day(ii),TZ_events.hour(ii),TZ_events.min(ii),floor(TZ_events.secs(ii)));
        tmp_lst = dir([DFolder,tmp_filename,'*']);
        tmp_filenames = {tmp_lst.name}';
        tmp_filedir = {tmp_lst.folder}';
        new_filenames = [new_filenames;tmp_filenames];
        new_filedir = [new_filedir;tmp_filedir];
    end
    
    listWithFolder = fullfile(new_filedir,new_filenames);
    listNoFolder = new_filenames;
    
end