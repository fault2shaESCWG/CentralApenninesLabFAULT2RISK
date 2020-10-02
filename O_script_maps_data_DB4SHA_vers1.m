% this code maps data stored in the Fault2SHA_CentralApennines_Database.xls
% and coordinates of MasterFaults according to the files given in the
% folder MasterFaults_lonlat
% and produces a map .......


% The user defines:
% 1) the maximum distance for which a point with slip rate is associated to the fault
% 2) the prefix of the output files name

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% start
clear all
clc
close all
warning('off','all')
addpath ('INPUT/','INPUT/MasterFaults_lonlat/')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% make output directory
pathout1 = fullfile('WORKING_DIRECTORY_A1B1C1_10km','Visualization','figure');

 
if isdir(pathout1)==0
mkdir (pathout1)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% USER OPTIONS
maxdiffUTM=100; % in meters, specify the max difference between two vertexes when resampling the masterfault trace
dmax = 500; % in meters, specify the maximum distance to associate a point to a masterfault
% colors for figures
limitisliprate = [0.0 0.1 0.5 1.0 3.0];
coloreslip = [.5 .5 .5; 0 1 0; 0 0 1; 1 0 0];
% fault certainActivity scale
values = 1:4;
coloreactivity = [1 0 0; 0 0 1; 0 1 0; .5 .5 .5];
latlim=([41.6 43.1]);
lonlim=([12.7 14.3]);


%% make figure
figure(1)
ax = usamap(latlim,lonlim);
setm(ax,'MapProjection','mercator', 'MapLatLimit',latlim,'MapLonLimit',lonlim,'FlineWidth',0.7,'FontSize',6);

hold on
%title('MasterFault traces and slip rate data')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% read the DB-excel format
[fault_data,masterfault_names,~] = xlsread ('Fault2SHA_CentralApennines_Database.xlsx','Fault');
[masterfault_data,masterfault_names2,~] = xlsread ('Fault2SHA_CentralApennines_Database.xlsx','MasterFault');

[sliprate_data,sliprate_TXT5,~] = xlsread ('Fault2SHA_CentralApennines_Database.xlsx','SlipRate');
[localgeomKin_data,localgeomKin_TXT6,~] = xlsread ('Fault2SHA_CentralApennines_Database.xlsx','LocalGeometryKinematics');

% faults name listed in the DB
masterfaults_all = masterfault_names(2:end,4);
masterfaults = unique(masterfaults_all);

% number of masterfaults in the DB
n_masterfault = size(masterfaults,1);
fprintf('you have %i masterfaults in the DB\n', n_masterfault)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

% position of the fault in the excel-Sheets of DB 
h1 = find(strcmp(masterfault_names2(:,2),masterfault))-1; % position-header
h5 = find(strcmp(sliprate_TXT5(:,4),masterfault))-1; % position-header
h6 = find(strcmp(localgeomKin_TXT6(:,4),masterfault))-1; % position-header

faultActivity = masterfault_data(h1,7);


% calculate the average dip from LocalGeometryKinematics 
average_dip = round(nanmean(localgeomKin_data(h6,8)),0); % average of dip along the entire FAULT, not for section
average_dip(isnan(average_dip))=60;


% CHECK if the order of the vertexes of the fault are given according to right-hand rule

strike = masterfault_data(h1,3);
az = azimuth(coordinateWGS(1,2),coordinateWGS(1,1),coordinateWGS(end,2),coordinateWGS(end,1));

d1 = az - strike;
d2 = az - (strike+180);
if abs(d2) < abs(d1)
    coordinateWGS = flipud(coordinateWGS);
end
[coordinateUTM(:,1),coordinateUTM(:,2),utmzonekml] = deg2utm(coordinateWGS(:,2),coordinateWGS(:,1));


% extract slip rate coordinates from the DB
period = sliprate_data(h5,[36]);
slipratecoordinateWGS = sliprate_data(h5,[16,17]);
slipratecoordinateUTM = sliprate_data(h5,[13,14]);
% calculate slip rate from slip and Throw in the DB
slipratePreferred     = sliprate_data(h5,24)./period; %
sliprateError         = sliprate_data(h5,[25,26])./[period, period];
SlipFromThrowPreferred        = (sliprate_data(h5,27)./sind(average_dip))./period;
SlipFromThrowError            = (sliprate_data(h5,[28,29])./sind(average_dip))./[period, period];

% if there are NaN values of slip rates in the database, (example
% pizzalto,MtVettore)
% we use throw and average dip to calculate the slip

hnan=[];
hnan = find(isnan(slipratePreferred));
slipratePreferred(hnan,:) = SlipFromThrowPreferred(hnan,:);
sliprateError(hnan,:)=SlipFromThrowError(hnan,:);
hnan=[];
% if you still have NaN then we remove that point!
hnan = find(isnan(slipratePreferred));
if ~isempty(hnan)
slipratecoordinateWGS(hnan,:)=[];
slipratecoordinateUTM(hnan,:) = [];
slipratePreferred(hnan,:) = [];
sliprateError(hnan,:)=[];
end



if sum(slipratePreferred >0) % at least a value of slip rate
    
% resample the fault trace
% the resampled trace is used to attribute point with slip rate to the
% fault

[y_fault_resUTM, x_fault_resUTM] = interpm(coordinateUTM(:,2),coordinateUTM(:,1),maxdiffUTM);
%%% check to remove double points
check_x=diff(x_fault_resUTM);
check_y=diff(y_fault_resUTM);

repeated_points=find(check_x<=1 & check_y<=1); % remove points with <1m distance to each other
x_fault_resUTM(repeated_points)=[];
y_fault_resUTM(repeated_points)=[]; 


% cumulative resampled-fault trace length
d_trace_resUTM =[];
 for i = 1:(length(x_fault_resUTM)-1)
     d_trace_resUTM(i+1) = sqrt((x_fault_resUTM(i,1)-x_fault_resUTM((i+1),1))^2+...
            (y_fault_resUTM(i,1)-y_fault_resUTM((i+1),1))^2);
 
 end
 
 resfault_cumsum_length =[];
 resfault_cumsum_length = cumsum(d_trace_resUTM);
 resfault_length = resfault_cumsum_length(end);


 

% associate the nearest point of the resampled fault trace
    d=[];min_d=[];
    for i = 1:size(slipratePreferred,1)
        temp_dist =[];
        temp_dist = sqrt((slipratecoordinateUTM(i,1)-x_fault_resUTM).^2+...
                         (slipratecoordinateUTM(i,2)-y_fault_resUTM).^2 );
                     
        d(i,1) = find(temp_dist==min(temp_dist),1,'first'); %save the ordinal position 
        min_d(i,1) =min(temp_dist);
    end

   % remove points with distance >min_d(m)
       d_outlimit = find(min_d>dmax);
       slipratePreferred(d_outlimit) = [];
       sliprateError(d_outlimit,:) = [];
       slipratecoordinateWGS(d_outlimit,:) = [];
       slipratecoordinateUTM(d_outlimit,:) = [];
       d(d_outlimit) = []; 


   if sum(slipratePreferred >0) % again, after removing points check if at least a positive slip rate exists   
    
distance_sliprate = resfault_cumsum_length(d);
% check if two or more meausures are on the same point
u = unique(distance_sliprate);
if length(u) < length(distance_sliprate)
    fprintf([masterfault,' there are 2 or more measures at the same location (used the mean) >>\n'])
    distanza_sliprate_2   =[];slipratePreferred_2   =[];sliprateError_2       =[];
    slipratecoordinateWGS_2 =[];slipratecoordinateUTM_2 =[];   
    for iu = 1:length(u)
       hu = find( distance_sliprate == u(iu));
       distanza_sliprate_2(iu) = mean(distance_sliprate(hu));
       slipratePreferred_2(iu,1) = mean(slipratePreferred(hu));
       sliprateError_2(iu,:) = mean(sliprateError(hu,:),1);
       slipratecoordinateWGS_2(iu,:) = mean(slipratecoordinateWGS(hu,:),1);
       slipratecoordinateUTM_2(iu,:) = mean(slipratecoordinateUTM(hu,:),1);
    end
    distance_sliprate       =[];slipratePreferred   = [];sliprateError       = [];
    slipratecoordinateWGS   = [];slipratecoordinateUTM = [];
    
    distance_sliprate       = distanza_sliprate_2;
    slipratePreferred       = slipratePreferred_2;
    sliprateError           = sliprateError_2;
    slipratecoordinateWGS   = slipratecoordinateWGS_2;
    slipratecoordinateUTM   = slipratecoordinateUTM_2;
    
end
   end
end

   %simplify section trace with a tollerance of tol degree about
    %example 0.005 is about 500meters
    
    [lat1, lon1, cerr, tol] = reducem(coordinateUTM(:,1),coordinateUTM(:,2),500);
    simplified_trace = [];
    simplified_trace = [lat1,lon1];
    
    utmzone = repmat(utmzonekml(1,:),size(simplified_trace,1),1);
    [Lat,Lon] = utm2deg(simplified_trace(:,1),simplified_trace(:,2),utmzone);
    
    Lat_simplifiedtrace = [Lat_simplifiedtrace;Lat;NaN];
    Lon_simplifiedtrace = [Lon_simplifiedtrace;Lon;NaN];
 

plotm(Lat_simplifiedtrace,Lon_simplifiedtrace,'-','color',coloreactivity(faultActivity,:))

if ~isempty(slipratePreferred)

[srcount,srl,srpos]=histcounts(slipratePreferred,limitisliprate);
% for i = 1:size(slipratecoordinateWGS,1)
% pmarker= plotm (slipratecoordinateWGS(i,2),slipratecoordinateWGS(i,1),'o','MarkerSize',4,...
%     'MarkerEdgeColor','k','MarkerFaceColor',(coloreslip(srpos(i),:)));
% 
% 
% end
hs=scatterm(slipratecoordinateWGS(:,2),slipratecoordinateWGS(:,1),12,slipratePreferred(:,1),'filled');

 hs.Children.MarkerFaceAlpha = .5;
 hs.Children.MarkerEdgeColor = 'k';
 colormap(coloreslip)
end


    end
end
% legend ACTIVITY
yl = 42.95;xl = 13.4;xl1 =13.50; xl2 =13.55;xl3 = 13.56;
textm(yl,xl,'Activity Scale','FontSize',6,'FontAngle','Italic')
for as = values
plotm([(yl-as*0.06),(yl-as*0.06)],[xl1,xl2],'-','color',coloreactivity(as,:));
textm((yl-as*0.06),xl3,num2str(as),'FontSize',6)
end

% legend SLIPRATES
yl = 42.95;xl = 13.8;xl1 =13.85; xl2 =13.9;
textm(yl,xl,'Slip rate (mm/a)','FontSize',6,'FontAngle','Italic')
labelsliprate = {'<0.1','0.1-0.5','0.6-1.0','>1.0'};
for srlim = 1:size(coloreslip,1)
plotm ((yl-srlim*0.06),xl1,'o','MarkerSize',3,...
    'MarkerEdgeColor','k','MarkerFaceColor',(coloreslip(srlim,:)))
textm((yl-srlim*0.06),xl2,labelsliprate(srlim),'FontSize',6)
end

% inset
h2 = axes('pos',[.31 .24 .2 .15]);
h2 = worldmap([35 46],[5 21]);
land = shaperead('landareas.shp', 'UseGeoCoords', true);
geoshow([land.Lat],[land.Lon])
ppatchm = patchm([latlim(1);latlim(1);latlim(2);latlim(2)],[lonlim(1);lonlim(2);lonlim(2);lonlim(1)],1)
ppatchm.FaceColor= 'r';
setm(h2, 'FFaceColor','w','FlineWidth',0.7)
mlabel; plabel; gridm % toggle off
hold off


saveas(1,fullfile(pathout1,'Maps_DB4SHA.png'),'tiff')
