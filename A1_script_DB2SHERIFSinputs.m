% this code reads input data stored in the
% Fault2SHA_CentralApennines_Database.xls % (https://doi.pangaea.de/10.1594/PANGAEA.922582)
% and coordinates of MainFaults according to the files given in the folder MainFaults_lonlat
% and produces SHERIFS input files.

% in the USER OPTIONS sections:

% The user has 2 options to calculate the average slip rate of a fault:
% a) arithmetic mean of the slip rate collected along the fault
% b) integral average assuming to have a slip rate = 0 at the tips of the fault

% The user defines:
% 1) the maximum distance for which a point with slip rate is associated to the fault
% 2) the prefix of the output files name
% 3) other variables required by SHERIFS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this code was written by Visini F. and Scotti O.
% please cite as: Scotti et al. (2020). 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
clc
close all
warning('off','all')
addpath ('INPUT/','INPUT/MainFaults_lonlat/')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% USER OPTIONS
extrapolate_slip_rate_option =1; % 1=tip to 0mm/year, 2 = no zero
modelname = 'ModelMultiFaultsTest'
maxdiffUTM= 0.1; % in km, specify the max difference between two vertexes when resampling the Mainfault trace
dmax = 0.5; % in km, specify the maximum distance to associate a point to a Mainfault
sections_length_input = 10; %km
maxdip = 55; % maximum dip to be assigned to the Mainfaults
USD = 0; %km
LSD = 10; %km
kin = 'N'; 
Domains = 'Active_Shallow_Crust';
Shmod = 30; % shear modulus
% colors for figures
limitisliprate = [0.0 0.1 0.5 1 3];
coloreslip = [.5 .5 .5; 0 1 0; 0 0 1; 1 0 0];
% fault certainActivity scale
values = 1:4;
coloreactivity = [1 0 0; 0 0 1; 0 1 0; .5 .5 .5];
% make output directory
mainpath = 'WORKING_DIRECTORY_A1B1C1_10km'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pathout1 = fullfile(mainpath,'Visualization','figure');
pathout2 = fullfile('A_SHERIFS_CAD','data','CAD_optionA1B1C1');
pathout3 = fullfile(mainpath,'Visualization','Data4Maps');
 
if isdir(pathout1)==0
mkdir (pathout1)
end
if isdir(pathout2)==0
mkdir (pathout2)
end
if isdir(pathout3)==0
mkdir (pathout3)
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%set units in kilometers
maxdiffUTM= maxdiffUTM*1000 ; 
dmax = dmax*1000;
sections_length_input = sections_length_input*1000; 

% read the DB-excel format
[fault_data,Mainfault_names,~] = xlsread ('Fault2SHA_CentralApennines_Database.xlsx','Fault');
[Mainfault_data,Mainfault_names2,~] = xlsread ('Fault2SHA_CentralApennines_Database.xlsx','MainFault');
[sliprate_data,sliprate_TXT5,~] = xlsread ('Fault2SHA_CentralApennines_Database.xlsx','SlipRate');
[localgeomKin_data,localgeomKin_TXT6,~] = xlsread ('Fault2SHA_CentralApennines_Database.xlsx','LocalGeometryKinematics');

[~,Mainfault_selection,~] = xlsread ('Fault2SHA_CentralApennines_Database.xlsx','MainFaultSelection','E7:E12');
R = Mainfault_selection;
display ('You have selected the following Main fault options')
R(~cellfun('isempty',R))
% faults name listed in the DB
Mainfaults_all = Mainfault_names(2:end,4);
Mainfaults = unique(Mainfaults_all);

% number of Mainfaults in the DB
n_Mainfault = size(Mainfaults,1);
fprintf('you have %i Mainfaults in the DB\n', n_Mainfault)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% open output files 
fid_prop = fopen(fullfile(pathout2,'Faults_properties.txt'),'w');
fid_geom = fopen(fullfile(pathout2,'Faults_geometry.txt'),'w');
fid_momentrate = fopen(fullfile(pathout2,'Faults_momentrate.txt'),'w');
fid_usedsliprate = fopen(fullfile(pathout3,'USED_sliprateDP.txt'),'w');
fid_usedMainfaults = fopen(fullfile(pathout3,'USED_Mainfaults.txt'),'w');

formatprop = '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n'; 
formatgeom = '%s\t%s\t%s\t%s\t%s\n';
formatmoment = '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n';
formatusedsliprate = '%s %s %s %s %s\n';
formatusedMainfaults = '%s\t%s\n';

fprintf(fid_prop,'model_name fault_name dip oriented kinematics upper_seismo_depth lower_seismo_depth slip_rate_min slip_rate_moy slip_rate_max Domain shear_modulus\n');
fprintf(fid_geom,'model_name fault_name longitude latitude type_of_fault\n');
fprintf(fid_momentrate,'model_name fault_name dip oriented kinematics upper_seismo_depth lower_seismo_depth slip_rate_min slip_rate_moy slip_rate_max Domain shear_modulus sectionlength MRmin MRmean MRmax\n');
fprintf(fid_usedsliprate,formatusedsliprate,'lon', 'lat','preferred', 'min','max');
fprintf(fid_usedMainfaults,formatusedMainfaults,'name', 'activityscale');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialize variables
All_tip_of_sections =[];
IDsection =0;

for nf = 1:size(Mainfaults,1)% start loop for each Mainfault
    sections_id_files  = [];
    LatSections = [];
    LonSections = [];
    Mainfault =[];
    Mainfault = char(Mainfaults(nf,:));
    % if Mainfaults exists in the folder MastreFaults_lonlat load coordinates of the Mainfaults
    in_name_Mainfault = fullfile('MainFaults_lonlat',strcat(Mainfault,'.txt'));
    if exist(in_name_Mainfault) ~= 2 % if a file exist the value is 2
      fprintf(['there are no coordinates of the ', Mainfault,' in the folder\n'])
    else  % use this MainFault
   
        coordinateUTM=[];coordinateWGS=[];
        coordinateWGS = load(in_name_Mainfault);

%% position of the fault in the excel-Sheets of DB 
h1 = find(strcmp(Mainfault_names2(:,2),Mainfault))-1; % position-header
h5 = find(strcmp(sliprate_TXT5(:,4),Mainfault))-1; % position-header
h6 = find(strcmp(localgeomKin_TXT6(:,4),Mainfault))-1; % position-header

faultActivity = Mainfault_data(h1,7);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf(fid_usedMainfaults,formatusedMainfaults, Mainfault,num2str(faultActivity));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate the average dip from LocalGeometryKinematics 
average_dip = round(nanmean(localgeomKin_data(h6,8)),0); % average of dip along the entire FAULT, not for section
average_dip(isnan(average_dip))= 60;
% if a maximum dip is given as input then it is used
if exist('maxdip','var')
    average_dip(average_dip > maxdip) = maxdip;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CHECK if the order of the vertexes of the fault are given according to right-hand rule

strike = Mainfault_data(h1,3);
az = azimuth(coordinateWGS(1,2),coordinateWGS(1,1),coordinateWGS(end,2),coordinateWGS(end,1));

d1 = az - strike;
d2 = az - (strike+180);
if abs(d2) < abs(d1)
    coordinateWGS = flipud(coordinateWGS);
end
[coordinateUTM(:,1),coordinateUTM(:,2),utmzonekml] = deg2utm(coordinateWGS(:,2),coordinateWGS(:,1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set the labels of the dipping  N,E,S,W
if (strike >0 & strike <= 30) | (strike >330 & strike <= 360 )
oriented = 'E';
elseif strike >30 & strike <= 120 
oriented = 'S';
elseif strike >120 & strike <= 210 
oriented = 'W';
elseif strike >210 & strike <= 330 
oriented = 'N';
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% extract slip rate coordinates from the DB
period = sliprate_data(h5,[36]);
slipratecoordinateWGS = sliprate_data(h5,[16,17]);
slipratecoordinateUTM = sliprate_data(h5,[13,14]);
% calculate slip rate from slip and Throw in the DB
slipratePreferred     = sliprate_data(h5,24)./period; %
sliprateError         = sliprate_data(h5,[25,26])./[period, period];
SlipFromThrowPreferred        = (sliprate_data(h5,27)./sind(average_dip))./period;
SlipFromThrowError            = (sliprate_data(h5,[28,29])./sind(average_dip))./[period, period];

% if there are NaN values of slip rates in the database, (example pizzalto,MtVettore)
% we use throw and average dip to calculate the slip

hnan=[];
hnan = find(isnan(slipratePreferred));
slipratePreferred(hnan,:) = SlipFromThrowPreferred(hnan,:);
sliprateError(hnan,:)=SlipFromThrowError(hnan,:);
hnan=[];
%% if you still have NaN then we remove that point!
hnan = find(isnan(slipratePreferred));
if ~isempty(hnan)
slipratecoordinateWGS(hnan,:)=[];
slipratecoordinateUTM(hnan,:) = [];
slipratePreferred(hnan,:) = [];
sliprateError(hnan,:)=[];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if sum(slipratePreferred >0) % at least a value of slip rate
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% cumulative fault trace length
for i = 1:(size(coordinateUTM,1)-1)
    d_trace(i+1) = sqrt((coordinateUTM(i,1)-coordinateUTM((i+1),1))^2+...
           (coordinateUTM(i,2)-coordinateUTM((i+1),2))^2);

end
fault_cumsum_length = cumsum(d_trace);
fault_length = fault_cumsum_length(end);


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

% associate the point with a measure of slip rate to nearest vertex of the resampled Mainfault trace
    d=[];min_d=[];
    for i = 1:size(slipratePreferred,1)
        temp_dist =[];
        temp_dist = sqrt((slipratecoordinateUTM(i,1)-x_fault_resUTM).^2+...
                         (slipratecoordinateUTM(i,2)-y_fault_resUTM).^2 );
                     
        d(i,1) = find(temp_dist==min(temp_dist),1,'first'); %save the ordinal position 
        min_d(i,1) =min(temp_dist);
    end
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % remove points with distance > min_d(m)
       d_outlimit = find(min_d>dmax);
       slipratePreferred(d_outlimit) = [];
       sliprateError(d_outlimit,:) = [];
       slipratecoordinateWGS(d_outlimit,:) = [];
       slipratecoordinateUTM(d_outlimit,:) = [];
       d(d_outlimit) = []; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   if sum(slipratePreferred >0) % again, after removing points check if at least a positive slip rate exists   
    
distance_sliprate = resfault_cumsum_length(d);
% check if two or more meausures are on the same point
u = unique(distance_sliprate);
if length(u) < length(distance_sliprate)
    fprintf([Mainfault,' there are 2 or more measures at the same location (used the mean) >>\n'])
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% loop to write points
for l_sr = 1:size(slipratecoordinateWGS,1)
fprintf(fid_usedsliprate,formatusedsliprate,num2str(slipratecoordinateWGS(l_sr,1)),num2str(slipratecoordinateWGS(l_sr,2)),...
                            num2str(slipratePreferred(l_sr)),num2str(sliprateError(l_sr,1)),num2str(sliprateError(l_sr,2)));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate average slip rate for the fault

    if extrapolate_slip_rate_option == 1

        if (min(distance_sliprate) == resfault_cumsum_length(1)) & (max(distance_sliprate) ~= resfault_cumsum_length(end))
                slipratecoordinateWGS = [slipratecoordinateWGS;coordinateWGS(end,:)];
                distance_sliprate   = [distance_sliprate,resfault_cumsum_length(end)];
                slipratePreferred   = [slipratePreferred;0];
                sliprateError       = [sliprateError;0,0];
        elseif (min(distance_sliprate) ~=resfault_cumsum_length(1)) & (max(distance_sliprate) == resfault_cumsum_length(end))
                slipratecoordinateWGS = [coordinateWGS(1,:);slipratecoordinateWGS];
                distance_sliprate   = [resfault_cumsum_length(1),distance_sliprate];
                slipratePreferred   = [0;slipratePreferred];
                sliprateError       = [0,0;sliprateError];    
        else   
            
    slipratecoordinateWGS = [coordinateWGS(1,:);slipratecoordinateWGS;coordinateWGS(end,:)];
    distance_sliprate   = [resfault_cumsum_length(1),distance_sliprate,resfault_cumsum_length(end)];
    slipratePreferred   = [0;slipratePreferred;0];
    sliprateError       = [0,0;sliprateError;0,0];
        end
        
% interpolate slip rate to build a profile        
    slipRateProfilePreferred = interp1(distance_sliprate,slipratePreferred,resfault_cumsum_length);
    slipRateProfileMin = interp1(distance_sliprate,sliprateError(:,1),resfault_cumsum_length);
    slipRateProfileMax = interp1(distance_sliprate,sliprateError(:,2),resfault_cumsum_length);
% integral average of the slip rate profile 
    a= min(distance_sliprate);
    b= max(distance_sliprate);
    AverageSlipRatePreferred=(1/(b-a))*trapz(distance_sliprate,slipratePreferred);
    AverageSlipRateMin=(1/(b-a))*trapz(distance_sliprate,sliprateError(:,1));
    AverageSlipRateMax=(1/(b-a))*trapz(distance_sliprate,sliprateError(:,2));
 
    elseif extrapolate_slip_rate_option == 2
 % mean of the data, no interpolation
    
    AverageSlipRatePreferred = mean(slipratePreferred);
    AverageSlipRateMin = mean(sliprateError(:,1));
    AverageSlipRateMax = mean(sliprateError(:,2));
    end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% define the sections
 
 n_section = round(resfault_length/sections_length_input);
 n_section (n_section ==0)=1;
 n_section (n_section ==2)=3;
 
 new_passo = resfault_length/n_section;
    limiti =[];
    limiti = new_passo:new_passo:resfault_length;
ind=[];
 for i = 1:length(limiti)
    ind(i) = find(abs( resfault_cumsum_length-limiti(i)) == min(abs( resfault_cumsum_length-limiti(i))));
end
ind = [1,ind];

% compute parameters of sections and write files 

    tip_of_sections =  [];% save tip of sections for plot
    meanslipratesectionPreferred = [];
    meanslipratesectionMin = [];
    meanslipratesectionMax = [];
    l_section = [];
for i =1:(length(ind)-1)
    section = [(x_fault_resUTM(ind(i):ind(i+1))), (y_fault_resUTM(ind(i):ind(i+1)))];
    
 l_section = [l_section; resfault_cumsum_length(ind(i+1))-resfault_cumsum_length(ind(i))];
 %%%%%%%%%%%%%%%%%%%%%
    %simplify section trace with a tolerance of 'tol' degree to reduce number of points defining a trace 
    %example tol=0.005 is roughly 500 meters
    
    [lat1, lon1, cerr, tol] = reducem(section(:,1), section(:,2),500);
    section = [];
    section = [lat1,lon1];
    
    utmzone = repmat(utmzonekml(1,:),size(section,1),1);
    [Lat,Lon] = utm2deg(section(:,1),section(:,2),utmzone);
    
    LatSections = [LatSections;Lat;NaN];
    LonSections = [LonSections;Lon;NaN];
  %%%%%%%%%%%%%%%%%%%%%%
       
    % save coordinates of the section in the fault_geom
    for j = 1: length(Lat)
    fprintf(fid_geom,formatgeom, modelname,strcat(Mainfault,'_',num2str(i)),num2str(Lon(j)), num2str(Lat(j)), 'sf');
    
    end 
    % save tip of sections for plot

    tip_of_sections =  [tip_of_sections;lon1(1) lat1(1); lon1(end), lat1(end)];
    All_tip_of_sections =[All_tip_of_sections; tip_of_sections];
    % save coordinate UTM of the sections
    IDsection = IDsection +1;
    sections_id_files = [sections_id_files;IDsection];
     UTMsectionsLON(1:size(section,1),IDsection)=section(:,1);
     UTMsectionsLAT(1:size(section,1),IDsection)=section(:,2);

     if extrapolate_slip_rate_option == 1
a=[];b=[]; 
a = resfault_cumsum_length(ind(i));
b = resfault_cumsum_length(ind(i+1));
meanslipratesectionPreferred(i,1) = (1/(b-a))*trapz(resfault_cumsum_length(ind(i):ind(i+1)),slipRateProfilePreferred(ind(i):ind(i+1))); %mean value theorem
meanslipratesectionMin(i,1) = (1/(b-a))*trapz(resfault_cumsum_length(ind(i):ind(i+1)),slipRateProfileMin(ind(i):ind(i+1))); %mean value theorem
meanslipratesectionMax(i,1) = (1/(b-a))*trapz(resfault_cumsum_length(ind(i):ind(i+1)),slipRateProfileMax(ind(i):ind(i+1))); %mean value theorem
     elseif extrapolate_slip_rate_option == 2
a=[];b=[]; 
a = resfault_cumsum_length(ind(i));
b = resfault_cumsum_length(ind(i+1));
 slip_in_section = [];
 ind_slip_in_section = find(distance_sliprate >=a & distance_sliprate <=b);
 if ~isempty(ind_slip_in_section)
 meanslipratesectionPreferred(i,1) = mean(slipratePreferred(ind_slip_in_section,1));
 meanslipratesectionMin(i,1) = mean(sliprateError(ind_slip_in_section,1));
 meanslipratesectionMax(i,1) = mean(sliprateError(ind_slip_in_section,2));
 elseif isempty(ind_slip_in_section)
    meanslipratesectionPreferred(i,1) = NaN;
    meanslipratesectionMin(i,1) = NaN;
    meanslipratesectionMax(i,1) = NaN;
 end  
     end
 
end

   %%%% treatement of NaN depending on position in the trace
    i_nan = [];i_not_nan = [];
    i_nan = find(isnan(meanslipratesectionPreferred));
    i_not_nan = find(~isnan(meanslipratesectionPreferred));
 
    pos = [];
    
 for pos = 1:length(i_nan)
    if i_nan(pos) ==1
     meanslipratesectionPreferred(i_nan(pos)) =  meanslipratesectionPreferred(i_not_nan(1)) ;
     meanslipratesectionMin(i_nan(pos)) =  meanslipratesectionMin(i_not_nan(1)) ;
     meanslipratesectionMax(i_nan(pos)) =  meanslipratesectionMax(i_not_nan(1)) ;
     
    elseif i_nan(pos) == length(meanslipratesectionPreferred)
        meanslipratesectionPreferred(i_nan(pos)) =  meanslipratesectionPreferred(i_not_nan(end)) ;
        meanslipratesectionMin(i_nan(pos)) =  meanslipratesectionMin(i_not_nan(end)) ;
        meanslipratesectionMax(i_nan(pos)) =  meanslipratesectionMax(i_not_nan(end)) ;
    else
        
        sx = i_not_nan(find(i_not_nan < i_nan(pos),1,'last'));
        dx = i_not_nan(find(i_not_nan > i_nan(pos),1,'first'));
        meanslipratesectionPreferred(i_nan(pos)) = mean([meanslipratesectionPreferred(sx),meanslipratesectionPreferred(dx) ]);
        meanslipratesectionMin(i_nan(pos)) = mean([meanslipratesectionMin(sx),meanslipratesectionMin(dx) ]);
        meanslipratesectionMax(i_nan(pos)) = mean([meanslipratesectionMax(sx),meanslipratesectionMax(dx) ]);

    end
 end
 

 for i=1: (length(ind)-1)
  
   
   fprintf(fid_prop,formatprop, modelname,strcat(Mainfault,'_',num2str(i)),...
                            num2str(average_dip),oriented,kin,num2str(USD),num2str(LSD),...
                            num2str(meanslipratesectionMin(i,1),'%3.2f'),num2str(meanslipratesectionPreferred(i,1),'%3.2f'),num2str(meanslipratesectionMax(i,1),'%3.2f'),...
                            Domains,num2str(Shmod));
                        
   MR1 = l_section(i,1)*Shmod*1e9* (LSD-USD)/sind(average_dip)*1000*  meanslipratesectionMin(i,1)/1000;
   MR2 = l_section(i,1)*Shmod*1e9* (LSD-USD)/sind(average_dip)*1000*  meanslipratesectionPreferred(i,1)/1000;
   MR3 = l_section(i,1)*Shmod*1e9* (LSD-USD)/sind(average_dip)*1000*  meanslipratesectionMax(i,1)/1000;
   
   fprintf(fid_momentrate,formatmoment, modelname,strcat(Mainfault,'_',num2str(i)),...
                            num2str(average_dip),oriented,kin,num2str(USD),num2str(LSD),...
                            num2str(meanslipratesectionMin(i,1),'%3.2f'),num2str(meanslipratesectionPreferred(i,1),'%3.2f'),num2str(meanslipratesectionMax(i,1),'%3.2f'),...
                            Domains,num2str(Shmod),num2str(l_section(i,1)),num2str(MR1),num2str(MR2),num2str(MR3));                     

 
 end
% 
figure(1)
hold on
% map of the fault with data points of slip rate
set(gca,'Visible','off')
  
lonlim=[12.7 14.2];
latlim=[41.7 43.1];
fontsize=10;fontname='Arial';
axesm ('mercator','frame','on','FlineWidth',0.7,'MapLatLimit',latlim,'MapLonLimit',lonlim)


title('fault trace, sections and slip rate data')
plotm(coordinateWGS(:,2),coordinateWGS(:,1),'-','color',coloreactivity(faultActivity,:))
plotm(tip_of_sections(:,2),tip_of_sections(:,1),'s','MarkerSize',2,'MarkerFaceColor','none','MarkerEdgeColor','k')
sect = 0;

[srcount,srl,srpos]=histcounts(slipratePreferred,limitisliprate);
for i = 1:size(slipratecoordinateWGS,1)
pmarker= plotm (slipratecoordinateWGS(i,2),slipratecoordinateWGS(i,1),'o','MarkerSize',5,...
    'MarkerEdgeColor','k','MarkerFaceColor',(coloreslip(srpos(i),:)));
pmarker.Color(4) = 0.3;
end
plotm(LatSections,LonSections,'-','color',coloreactivity(faultActivity,:))

hold off
xlabel ('Latitude')
ylabel('Longitude')

figure(2)
leg=[];
hleg = 0;
hold on

if extrapolate_slip_rate_option ==1
    hleg = hleg+1;
leg(1,hleg) = plot(resfault_cumsum_length,slipRateProfilePreferred,'-','color',[.5 .5 .5],'LineWidth',2,'Display','Interpolated Preferred');
hleg = hleg+1;
leg(1,hleg) = plot(resfault_cumsum_length,slipRateProfileMin,':','color',[.5 .5 .5],'LineWidth',2,'Display','Interpolated Min Max');
                plot(resfault_cumsum_length,slipRateProfileMax,':','color',[.5 .5 .5],'LineWidth',2,'Display','InterpolatedMax');
hleg = hleg+1;
leg(1,hleg) =  line([distance_sliprate(1),distance_sliprate(end)],[AverageSlipRatePreferred,AverageSlipRatePreferred],'color','c','LineWidth',1,'Display','Avg Preferred');
hleg = hleg+1;
leg(1,hleg) =  plot([distance_sliprate(1),distance_sliprate(end)],[AverageSlipRateMin,AverageSlipRateMin],':','color','c','LineWidth',1,'Display','Avg Min Max');
            plot([distance_sliprate(1),distance_sliprate(end)],[AverageSlipRateMax,AverageSlipRateMax],':','color','c','LineWidth',1,'Display','Avg Max');

elseif extrapolate_slip_rate_option ==2
   hleg = hleg+1;
leg(1,hleg) = line([resfault_cumsum_length(1),resfault_cumsum_length(end)],[AverageSlipRatePreferred,AverageSlipRatePreferred],'color','c','LineWidth',1,'Display','Avg Preferred');
   hleg = hleg+1;
leg(1,hleg) = plot([resfault_cumsum_length(1),resfault_cumsum_length(end)],[AverageSlipRateMin,AverageSlipRateMin],':','color','c','LineWidth',1,'Display','Avg Min Max');
  hleg = hleg+1;
leg(1,hleg) =  plot([resfault_cumsum_length(1),resfault_cumsum_length(end)],[AverageSlipRateMax,AverageSlipRateMax],':','color','c','LineWidth',1,'Display','Avg Min Max');

end

hleg1 = hleg+1;
hleg2 = hleg+2;
for i = 1: (length(meanslipratesectionPreferred))
    
        leg(1,hleg1) =    line([resfault_cumsum_length(ind(i)) resfault_cumsum_length(ind(i+1))],[meanslipratesectionPreferred(i),meanslipratesectionPreferred(i)],'color','m','LineWidth',1,'Display','Avg Sect Pref');
  
        leg(1,hleg2) =    plot([resfault_cumsum_length(ind(i)) resfault_cumsum_length(ind(i+1))],[meanslipratesectionMin(i),meanslipratesectionMin(i)],':','color','m','LineWidth',1,'Display','Avg Sect Min Max');
                          plot([resfault_cumsum_length(ind(i)) resfault_cumsum_length(ind(i+1))],[meanslipratesectionMax(i),meanslipratesectionMax(i)],':','color','m','LineWidth',1,'Display','Avg Sect Max');

end

hleg = hleg2+1;
leg(1,hleg)=errorbar(distance_sliprate,slipratePreferred,(slipratePreferred-sliprateError(:,1)),(sliprateError(:,2)-slipratePreferred),...
    'Marker','s','MarkerFaceColor','w','LineStyle','none', 'MarkerEdgeColor','k','Color','k','Display','DataPoints');
xlabel('distance along strike (km)')
ylabel ('mm/yr')
legend(leg,'location','best')

% add rough orientation of the fault strike
label1 = {'S','SW', 'NW', 'N','N','NW', 'SW', 'S'};
label2 = {'N','NE', 'SE', 'S','S','SE', 'NE', 'N'};
limits = [0 45 90 135 180 225 270 315 360];
label_pos = find(histcounts(strike,limits)==1);
text(-1500,AverageSlipRatePreferred,label1(label_pos));
text((resfault_length+500),AverageSlipRatePreferred,label2(label_pos));

xlim([-2000,resfault_length+1000])
ylim([0 3])

set(gca, 'XTick',[-2000:2000:(resfault_length+500)],'XTickLabel',{'',(0:2000:(resfault_length+500))/1000})


saveas(2,fullfile(pathout1,strcat(Mainfault,'_sliprates_profile_option',num2str(extrapolate_slip_rate_option),'_',date,'.png')),'png')
print(fullfile(pathout1,strcat(Mainfault,'_sliprates_profile_option',num2str(extrapolate_slip_rate_option),'_',date,'.tiff')),'-dtiff','-r600');

close(2);
   end
end
    end
end

saveas(1,fullfile(pathout1,strcat('map','_sliprates_profile_option',num2str(extrapolate_slip_rate_option),'_',date,'.png')),'tiff')
print(fullfile(pathout1,strcat('map','_sliprates_profile_option',num2str(extrapolate_slip_rate_option),'_',date,'.tiff')),'-dtiff','-r600');

fclose(fid_prop);
fclose(fid_geom);
fclose(fid_momentrate);
fclose(fid_usedMainfaults);
fclose(fid_usedsliprate);
