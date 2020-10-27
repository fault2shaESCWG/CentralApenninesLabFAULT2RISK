% map of section with mean slip rates from SHERPFS
clear all
clc
close all
warning('off','all')

limitisliprate = [0 0.1 0.5 1 3];
coloreslip = [.5 .5 .5; 0 1 0; 0 0 1; 1 0 0];
labelsliprate = {'<0.1','0.1-0.5','0.6-1.0','>1.0'};
latlim=([41.6 43.2]);
lonlim=([12.7 14.3]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%USER OPTIONS

mainpath = 'WORKING_DIRECTORY_A1B1C1_10km';
model_output = fullfile(mainpath,'Visualization');
sherifs_path1 = fullfile('A_SHERIFS_CAD','data','CAD_optionA1B1C1_10km');
sherifs_path2 = fullfile('A_SHERIFS_CAD','CAD_optionA1B1C1_10km','analysis','txt_files');

fault_prop = readtable(fullfile(sherifs_path1,'Faults_properties.txt'));
fault_geom = readtable(fullfile(sherifs_path1,'Faults_geometry.txt'));
fault_slip = readtable(fullfile(sherifs_path2,'mean_parameters_faults.txt'));
fault_nms = readtable(fullfile(sherifs_path2,'slip_rep_on_faults_all_data.txt'));


paleo = [55,11,24,26]; % CI VANNO I NUMERI DELLE SECTION ORA SONO NUMERI A CASO

 %% mappa section
 
 figure(1)

 hold on
 ax = usamap(latlim,lonlim);
 setm(ax,'MapProjection','mercator', 'MapLatLimit',latlim,'MapLonLimit',lonlim,'FlineWidth',0.7,'FontSize',6);
%

 
 for i = 1:size(fault_prop,1)
 % colored by INPUT slip rate
   slip_input = fault_slip.Var4(i);
   nms_output = fault_nms.Var22(i);
 slip_output = slip_input;
   %slip_output = slip_input - slip_input*nms_output/100;
   slip_pos = find(histcounts(slip_output,limitisliprate));  
id = fault_prop.Var2{i};
pos = find(strcmp(fault_geom.Var2,id)==1);
xf = fault_geom.Var3(pos);
yf = fault_geom.Var4(pos);


% PALEO
if sum(i == paleo)>=1
plotm(mean(yf),mean(xf),'Marker','pentagram','MarkerSize',10,'LineWidth',1,'MarkerFaceColor','none','MarkerEdgeColor',[.5 .5 .5])
end

plotm(yf,xf,'-','LineWidth',1,'color',coloreslip(slip_pos,:))

plotm(yf(1),xf(1),'s','MarkerSize',2,'LineWidth',1,'MarkerFaceColor','none','MarkerEdgeColor','k')
plotm(yf(end),xf(end),'s','MarkerSize',2,'LineWidth',1,'MarkerFaceColor','none','MarkerEdgeColor','k')


% numero section
   yc = mean(yf)+0.01;
   xc =mean(xf)+0.005;
   %testo1 = strcat('(',num2str(i),')');
   testo1 = strcat(num2str(i));
   textm(yc,xc,testo1,'FontSize',6,'BackgroundColor','none','Margin',0.001)


 end

% legend tip of sections
yl = 42.95;xl1 = 13.44;xl2 =13.49; 
plotm(yl,xl1,'s','MarkerSize',2,'LineWidth',1,'MarkerFaceColor','none','MarkerEdgeColor','k')
textm(yl,xl2,'tip of sections','FontSize',6)

% legend PALEO
yl = 42.90;xl1 = 13.44;xl2 =13.49; 
plotm(yl,xl1,'Marker','pentagram','MarkerSize',10,'LineWidth',1,'MarkerFaceColor','none','MarkerEdgeColor',[.5 .5 .5])
textm(yl,xl2,'section with paleoseimological data','FontSize',6)

% legend of SLIPRATES

legend1 = []
srleg = []
for srlim = 1:size(coloreslip,1)
srleg(srlim,1) = plotm (90,90,'-','LineWidth',1,'color',coloreslip(srlim,:),...
    'Display',char(labelsliprate(srlim)));
end
legend1 = legend(srleg)
legend1.Position = [0.38 0.25 0.07 0.12];
legend1.Box = 'off';
legend1.FontSize = 6;
title(legend1,'slip rate (mm/yr)')
legend1.Title.Visible = 'on';




hold off % termina la figura 1 mappa
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

saveas(1,fullfile(model_output,'figure',strcat('map_sections','_',date,'.tiff')),'tiff')

