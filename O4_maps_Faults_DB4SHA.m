% this code maps faults stored in the Fault2SHA_CentralApennines_Database.xls

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% start
clear all
clc
close all
warning('off','all')
addpath ('INPUT/','INPUT/MainFaults_lonlat/')
addpath ('SHAPEFILES/')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% USER OPTIONS
% colors for figures
limitisliprate = [0.0 0.1 0.5 1.0 3.0];
coloreslip = [.5 .5 .5; 0 1 0; 0 0 1; 1 0 0];
% fault certainActivity scale
values = 1:4;
coloreactivity = [255,0,0; 255,125,0; 255,195,0; 255,255,0]/255;
list = {'GioiaVecchio', 'Magnola','OvindoliPezza','SanBenedettoDeiMarsi',...
    'CastelDiIeri','Frattura','LagoDiScanno','Bazzano','Collebrincioni',...
    'Fossa','MtRuzza','MtStabiata','Paganica','SanDemetrioNeVestini',...
    'SanDemetrioEast','SanPioDelleCamere'};

latlim=([41.6 43.1]);
lonlim=([12.7 14.3]);
% make output directory
pathout1 = fullfile('WORKING_DIRECTORY_A1B1C1_10km','Visualization','figure');
%pathout1 = fullfile('WORKING_DIRECTORY_A2B2C2_10km','Visualization','figure');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isdir(pathout1)==0
mkdir (pathout1)
end
%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% read the traces from the shapefiles
input_faults = shaperead('faults.shp');
faultName = {input_faults.Name}';

% read the DB-excel format
db = 'Fault2SHA_CentralApennines_Database_2021_v1_xls2013.xlsx';

[fault_data,Mainfault_names,~] = xlsread (db,'Fault');
[Mainfault_data,Mainfault_names2,~] = xlsread (db,'MainFault');
FaultActivity = fault_data(:,8);
[sliprate_data,sliprate_TXT5,~] = xlsread (db,'SlipRate');
[localgeomKin_data,localgeomKin_TXT6,~] = xlsread (db,'LocalGeometryKinematics');

% number of traces in the DB
n_traces = size(faultName,1);
fprintf('you have %i faults in the DB\n', n_traces)


% extract slip rate coordinates from the DB
period = sliprate_data(:,[36]);
slipratecoordinateWGS = sliprate_data(:,[16,17]);
% calculate slip rate from slip and Throw in the DB
slipratePreferred     = sliprate_data(:,24)./period; %
sliprateError         = sliprate_data(:,[25,26])./[period, period];
ThrowratePreferred        = (sliprate_data(:,27))./period;
ThrowrateError            = (sliprate_data(:,[28,29]))./[period, period];

% if there are NaN values of slip rates in the database, (example pizzalto,MtVettore) we use throw and average dip to calculate the slip

hnan=[];
hnan = find(isnan(slipratePreferred));
% calculate the average dip from LocalGeometryKinematics to be used for
% computing slip rate from throw rate where necessary
for inan = 1:length(hnan)

h6 = find(strcmp(localgeomKin_TXT6(:,4),sliprate_TXT5(hnan(inan)+1,4)))-1; % position-header
average_dip = round(nanmean(localgeomKin_data(h6,8)),0); % average of dip along the entire FAULT, not for section
average_dip(average_dip>60 | isnan(average_dip)) =55;
slipratePreferred(hnan(inan),:) = ThrowratePreferred(hnan(inan),:)./sind(average_dip);
sliprateError(hnan(inan),:) = ThrowrateError(hnan(inan),:)./sind(average_dip);
end
hnan=[];
% if you still have NaN then we remove that point!
hnan = find(isnan(slipratePreferred));
if ~isempty(hnan)
slipratecoordinateWGS(hnan,:)=[];
slipratePreferred(hnan,:) = [];
sliprateError(hnan,:)=[];
end


%% make figure
figure(1)
ax = usamap(latlim,lonlim);
setm(ax,'MapProjection','mercator', 'MapLatLimit',latlim,'MapLonLimit',lonlim,...
     'FlineWidth',0.7,'FontSize',8,'MLabelLocation',1,'PLabelLocation',1,'MLabelRound',0,'PLabelRound',0,'LabelFormat','none');

hold on

for i_trace = 1:n_traces
trace_name = faultName(i_trace,:)
if ismember(trace_name,list)==1
    plotm(input_faults(i_trace,:).Y,input_faults(i_trace,:).X,'-','LineWidth',1.4,'color','k')
end
plotm(input_faults(i_trace,:).Y,input_faults(i_trace,:).X,'-','LineWidth',1.,'color',coloreactivity(FaultActivity(i_trace),:))

end
% hs=scatterm(slipratecoordinateWGS(:,2),slipratecoordinateWGS(:,1),12,slipratePreferred(:,1),'filled');
% hs.Children.MarkerFaceAlpha = .5;
% hs.Children.MarkerEdgeColor = 'k';
% colormap(coloreslip)

 
% legend ACTIVITY
yl = 42.95;xl = 13.56;xl1 =13.60; xl2 =13.65;xl3 = 13.66;
textm(yl,xl,'Activity Scale','FontSize',6,'FontAngle','Italic')
for as = values
plotm([(yl-as*0.06),(yl-as*0.06)],[xl1,xl2],'-','LineWidth',1.5,'color',coloreactivity(as,:));
textm((yl-as*0.06),xl3,num2str(as),'FontSize',6)
end

% legend SLIPRATES
% yl = 42.95;xl = 13.8;xl1 =13.88; xl2 =13.92;
% textm(yl,xl,'Slip rate (mm/a)','FontSize',8,'FontAngle','Italic')
% labelsliprate = {'<0.1','0.1-0.5','0.6-1.0','>1.0'};
% for srlim = 1:size(coloreslip,1)
% plotm ((yl-srlim*0.06),xl1,'o','MarkerSize',3,'LineWidth',0.5,...
%     'MarkerEdgeColor','k','MarkerFaceColor',(coloreslip(srlim,:)))
% textm((yl-srlim*0.06),xl2,labelsliprate(srlim),'FontSize',8)
% end

% inset
%h2 = axes('pos',[.31 .24 .2 .15]);
% h2 = axes('pos',[.32 .26 .18 .11]);
% h2 = worldmap([35 46],[5 21]);
% land = shaperead('landareas.shp', 'UseGeoCoords', true);
% geoshow([land.Lat],[land.Lon])
% ppatchm = patchm([latlim(1);latlim(1);latlim(2);latlim(2)],[lonlim(1);lonlim(2);lonlim(2);lonlim(1)],1);
% ppatchm.FaceColor= 'r';
% setm(h2, 'FFaceColor','w','FlineWidth',0.7)
% mlabel; plabel; gridm % toggle off
hold off

%% move label of meridians and parallels
set(findobj(ax.Children, 'Tag', 'MLabel'),'Units','points')          % convert label position from 'data' to 'points'
mlabels = findobj(ax.Children, 'Tag', 'MLabel');                     % find all labels
mlabelpos = get(findobj(ax.Children, 'Tag', 'MLabel'),'Position');    % get the positions of each label
for iL = 1 : length(mlabelpos)                                           % loop over each label
    mlabelpos{iL}(2) = mlabelpos{iL}(2) + 10;                             % add desired offset to the label position
    set(mlabels(iL),'Position',mlabelpos{iL})                            % set new label position
end
set(findobj(ax.Children, 'Tag', 'PLabel'),'Units','points')          % convert label position from 'data' to 'points'
plabels = findobj(ax.Children, 'Tag', 'PLabel');                     % find all labels
plabelpos = get(findobj(ax.Children, 'Tag', 'PLabel'),'Position');    % get the positions of each label
for iL = 1 : length(plabelpos)                                           % loop over each label
    plabelpos{iL}(1) = plabelpos{iL}(1) + 4;
    plabelpos{iL}(2) = plabelpos{iL}(2) + 1;                    % add desired offset to the label position
    set(plabels(iL),'Position',plabelpos{iL})                            % set new label position
end


saveas(1,fullfile(pathout1,'MapsFaults_DB4SHA.png'),'tiff')
print(fullfile(pathout1,'MapsFaults_DB4SHA.tiff'),'-dtiff','-r600');