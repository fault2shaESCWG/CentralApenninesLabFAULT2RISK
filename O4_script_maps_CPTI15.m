% this code maps data stored in the Fault2SHA_CentralApennines_Database.xls
% and coordinates of MasterFaults according to the files given in the
% folder MasterFaults_lonlat
% and CPTI15 .......


% The user defines:
% 1) the prefix of the output files name

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% start
clear all
clc
close all
warning('off','all')
addpath ('INPUT/','INPUT/MasterFaults_lonlat/')
addpath ('INPUT/','area/')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% make output directory
pathout1 = fullfile('WORKING_DIRECTORY_A1B1C1_10km','Visualization','figure');

 
if isdir(pathout1)==0
mkdir (pathout1)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% USER OPTIONS

latlim=([41.6 43.1]);
lonlim=([12.7 14.3]);

minimum_magnitude = 5;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% read CPTI15 and area
 cpti15 = readtable('CPTI15_extracted.csv');
 x_eq = cpti15.LonDef(cpti15.MwDef >= minimum_magnitude);
 y_eq = cpti15.LatDef(cpti15.MwDef >= minimum_magnitude);
 m_eq = cpti15.MwDef(cpti15.MwDef >= minimum_magnitude);

 area = shaperead(fullfile('INPUT','area','background.shp'));
 x_area = [area.X]';y_area = [area.Y]';
%% read the DB-excel format
[fault_data,masterfault_names,~] = xlsread ('Fault2SHA_CentralApennines_Database.xlsx','Fault');
[masterfault_data,masterfault_names2,~] = xlsread ('Fault2SHA_CentralApennines_Database.xlsx','MasterFault');

% faults name listed in the DB
masterfaults_all = masterfault_names(2:end,4);
masterfaults = unique(masterfaults_all);

% number of masterfaults in the DB
n_masterfault = size(masterfaults,1);
fprintf('you have %i masterfaults in the DB\n', n_masterfault)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% make figure
figure(1)
ax = usamap(latlim,lonlim);
setm(ax,'MapProjection','mercator', 'MapLatLimit',latlim,'MapLonLimit',lonlim,'FlineWidth',0.7,'FontSize',6);

hold on
%title('MasterFault traces and slip rate data')

%% loop for masterfaults
IDsection =0;

for nf = 1:size(masterfaults,1);
    Lat_simplifiedtrace = [];
    Lon_simplifiedtrace = [];
    masterfault =[];
    masterfault = char(masterfaults(nf,:));
    % if masterfaults exists in the folder MastreFaults_lonlat load coordinates of the masterfaults
    in_name_masterfault = fullfile('MasterFaults_lonlat',strcat(masterfault,'.txt'));
    if exist(in_name_masterfault) ~= 2 % if a file exist the value is 2
      fprintf(['there are no coordinates of the ', masterfault,' in the folder\n'])
    else  % use this MasterFault
        coordinateUTM=[];coordinateWGS=[];
        coordinateWGS = load(in_name_masterfault);

[coordinateUTM(:,1),coordinateUTM(:,2),utmzonekml] = deg2utm(coordinateWGS(:,2),coordinateWGS(:,1));


   %simplify the trace with a tollerance of tol degree about
    %example 0.005 is about 500meters
    
    [lat1, lon1, cerr, tol] = reducem(coordinateUTM(:,1),coordinateUTM(:,2),500);
    simplified_trace = [];
    simplified_trace = [lat1,lon1];
    
    utmzone = repmat(utmzonekml(1,:),size(simplified_trace,1),1);
    [Lat,Lon] = utm2deg(simplified_trace(:,1),simplified_trace(:,2),utmzone);
    
    Lat_simplifiedtrace = [Lat_simplifiedtrace;Lat;NaN];
    Lon_simplifiedtrace = [Lon_simplifiedtrace;Lon;NaN];
 

plotm(Lat_simplifiedtrace,Lon_simplifiedtrace,'-','LineWidth',1.5,'color',[0 0 0])

    end
end

%% add cpti15 and area and legend CPTI
parea = plotm(y_area,x_area);
y_eq_label = 43.0; x_eq_label = 13.5;
mc = 5:0.5:7.5;
for lmc = 1:(length(mc)-1)
    y_eq_label=y_eq_label-0.07;
    hmc = find(m_eq>=mc(lmc) & m_eq <mc(lmc+1));
    scatterm(y_eq(hmc),x_eq(hmc),mc(lmc).^2.5,'o');
    scatterm(y_eq_label,x_eq_label,mc(lmc)^2.5,'o');...
    textm(y_eq_label,x_eq_label+0.04,strcat(num2str(mc(lmc),'%.2f'),'<=Mw<',num2str(mc(lmc+1),'%.2f')),'fontsize',6);

end


saveas(1,fullfile(pathout1,'Maps_CPTI15.png'),'tiff')
print(fullfile(pathout1,'Maps_CPTI15'),'-depsc','-r600');