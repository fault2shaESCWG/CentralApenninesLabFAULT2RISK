% this code calculates all the distances betwenn all the sections in the fault system and 
% creates rupture scenarios for input to SHERIFS according to the user defined criteria
% Coordinates of the sections are produced by A1_scipt_DB2SHERIFSinput.m 

% In the USER OPTIONS section the user needs to define the maximum possible 
% distance between sections capable of rupturing in a single scenario 
% the maximum combination of sections making up a rupture
% Set the flag to 1 to visualze rupture scenarios 

clear all
clc
close all
warning('off','all')

E = wgs84Ellipsoid('meters');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% USER OPTIONS
% maximum possible distance to build a scenario of 2 sections is 5 km
maxD = 5;
% the maximum combination of sections making up a rupture
possible_max_combinations = 7;
% make figures
makefigure =0;
%define set ID
setID = 'Set_1';
% make output directory
mainpath = 'WORKING_DIRECTORY_A1B1C1_10km'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pathout1 = fullfile(mainpath,'Visualization','scenarios');
pathout2 = fullfile('A_SHERIFS_CAD','input','CAD_optionA1B1C1');
if isdir(pathout1)==0
mkdir (pathout1)
end
if isdir(pathout2)==0
mkdir (pathout2)
end

% output for SHERIFS
fid_scenari = fopen(fullfile(pathout2,'ruptures.txt'),'w');
fprintf(fid_scenari, ['set ', setID,'\n']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% READ INPUTS
pathin1 = fullfile('A_SHERIFS_CAD','data','CAD_optionA1B1C1');
fault_geom = readtable(fullfile(pathin1,'Faults_geometry.txt'));
section_name = unique(fault_geom.Var2);

IDsection = 1:size(section_name);

UTMsectionsLON(1:1000, 1:IDsection(end))=NaN;
UTMsectionsLAT(1:1000, 1:IDsection(end))=NaN;

for i=1:IDsection(end) 
    section = section_name(i,:);
    pos = find(strcmp(fault_geom.Var2,section));
    UTMsectionsLON(1:length(pos),i) = fault_geom.Var3(pos);
    UTMsectionsLAT(1:length(pos),i) = fault_geom.Var4(pos);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n_sections = IDsection(end);

coppie_all(1:n_sections,1:n_sections) = repmat(1:n_sections,n_sections,1);
 
coppie = triu(coppie_all);
A=0;
% compute the minimum distance between all possible couples of sections and
% create the matrix 'mind'
for i =1:n_sections
    A=A+1;
   xA =[];yA=[];
   xA= UTMsectionsLON(:,A);
   yA= UTMsectionsLAT(:,A);
   xA(isnan(xA)) = [];
   yA(isnan(yA)) = [];
    wh_sect = [];
    
    wh_sect =  coppie(i,:);
    wh_sect = wh_sect(wh_sect>0);
  
 for j = 1:  length(wh_sect);
  templong =   UTMsectionsLON(:,wh_sect(j));
  templong(isnan(templong))=[];
  templat =   UTMsectionsLAT(:,wh_sect(j));
  templat(isnan(templat))=[];
  
     d=[];
     for z = 1 :length(xA)
        [ARCLEN, AZ] = distance(xA(z),yA(z),templong,templat,E);
        d(z,1) = min(ARCLEN);
     end
   
   mind(i,wh_sect(j)) = min(d)/1000;% km 
   mind(wh_sect(j),i) = min(d)/1000;% km 
 end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% extract couples of sections according to user option criteria

z=0;
for i = 1 : size(mind,1)
    for j = i : size(mind,2)
        if (i~=j) & (mind(i,j) <= maxD)
            z = z+1;
    temp_scenari1(z,1:2) = [i,j];
        end
    end
end

% remove scenario obtained by combination: 1+2 2+1 etc etc
temp_scenari1=sort(temp_scenari1,2);
temp_scenari2 = unique(temp_scenari1,'rows');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% add the first set of scenario to the final matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
scenari(1:2000,1:possible_max_combinations)=NaN;
scenari(1:size(temp_scenari2,1),1:size(temp_scenari2,2))=temp_scenari2;
last_pos = find(isnan(scenari(:,1)),1,'first');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ini = 1;
for q = 1:(possible_max_combinations-2) 
    
%%%% combine three or more sections
temp_scenari1 =[];temp_scenari2 = [];temp_scenari3=[];
last_col =[];
last_col = find(~isnan(scenari((last_pos-1),:)),1,'last');
temp_scenari1 = scenari(ini:(last_pos-1),1:last_col);
%pause

pos = 1;
for i = 1 : size(temp_scenari1,1)
    for j = 1: n_sections
        if prod(j-temp_scenari1(i,:))~=0 % 0 distances are not counted
            % find the minimum distance of the matrix ' mind'
            dist =[];
           
            for k = 1:length(temp_scenari1(i,:))
             %    [temp_scenari1(i,:),j]
            dist(k) = mind(temp_scenari1(i,k),j);
            %pause
            end
            if min(dist) <= maxD
            temp_scenari = [temp_scenari1(i,:),j];
            temp_scenari2(pos,1:length(temp_scenari))=temp_scenari;
            pos = pos+1;
            end
        end
    end
end

% remove scenario obtained by combination: 1+2+3 2+1+3 etc etc
temp_scenari2=sort(temp_scenari2,2);
temp_scenari3=unique(temp_scenari2,'rows');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
scenari(last_pos:(last_pos-1+size(temp_scenari3,1)),...
        1:size(temp_scenari3,2))=temp_scenari3;
    ini = last_pos;
last_pos = find(isnan(scenari(:,1)),1,'first');
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
end

 if makefigure ==1
for i =1 :size(scenari,1)
    figure(2)
    hold on
    last_col = find(~isnan(scenari(i,:)),1,'last');
    for j = 1:last_col
    plot(UTMsectionsLON(:,scenari(i,j)), UTMsectionsLAT(:,scenari(i,j)))
    end
    hold off
    xlim([min(min(UTMsectionsLON)) max(max(UTMsectionsLON))]);
    ylim([min(min(UTMsectionsLAT)) max(max(UTMsectionsLAT))]);
    saveas(2,fullfile(pathout1,strcat(num2str(i),'_scenario_',date,'.png')),'png')
    close(2)
end
 end
 
 last_pos_save = find(~isnan(scenari(:,1)),1,'last');
 scenari = scenari(1:last_pos_save,:);
 
 save(fullfile(pathout1,'scenari.txt'),'scenari','-ascii')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% WRITE SCENARIOS 
last_pos_save = find(~isnan(scenari(:,1)),1,'last');
scenari = scenari(1:last_pos_save,:);

fprintf(['Number of rupture scenarios : ',num2str(size(scenari,1)), ' in ', setID,'\n'])
fprintf(['Maximum distance between sections : ',num2str(maxD),' \n']);
fprintf(['Input max combinations: ',num2str(possible_max_combinations),' \n']);
lscenari = sum(scenari(end,:) >0);
fprintf(['Output max combinations: ',num2str(lscenari),' \n']);
single_sections = unique(IDsection);
single_sections (isnan(single_sections))=[];
        
for i = 1: size(scenari,1)
    last_col =[];

    last_col = find((scenari(i,:))>0,1,'last');
    for j = 1:last_col
        row_section = find(IDsection == scenari(i,j));
        if j < last_col
        fprintf(fid_scenari,'%s ', strcat(section_name{row_section,:}));
        elseif j == last_col
        fprintf(fid_scenari, strcat(section_name{row_section,:},'\n'));    
            end
        end
    end