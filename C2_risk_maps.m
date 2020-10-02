% This script builds a map of fault sections colored coded according to the
% mean slip rates used in SHERIFS and risk estimates at all the
% localities of the study regions

clear all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% USER OPTIONS
fragility ='Rostietal.2020-L-type';%set fragility
taxonomy = 'pre1919';
%taxonomy = 'post1981';

mainpath = 'WORKING_DIRECTORY_A1B1C1_10km';
openquakepath = fullfile(mainpath,'OQoutputs');
model_output = fullfile(mainpath,'Visualization');
sherifs_path1 = fullfile('A_SHERIFS_CAD','data','CAD_optionA1B1C1');
sherifs_path2 = fullfile('A_SHERIFS_CAD','CAD_optionA1B1C1','analysis','txt_files');

OQ_RUN_ID = '3';% Number of Openquake run ID
fprintf(['Warning: You are using OQ_RUN_ID ',num2str(OQ_RUN_ID)]);

limitisliprate = [0 0.1 0.5 1 3];
coloreslip = [.5 .5 .5; 0 1 0; 0 0 1; 1 0 0];
labelsliprate = {'<0.1','0.1-0.5','0.6-1.0','>1.0'};
mycolors = [.5 .5 .5; 0 1 0; 0 0 1; 1 0 0];
latlim=([41.6 43.2]);
lonlim=([12.7 14.3]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

realization = readtable(fullfile(openquakepath,strcat('realizations_',num2str(OQ_RUN_ID),'.csv')));
nrlz= size(realization,1)-1; % number of scenarios A1B1C1, note that the branchID goes from 0 to n-1

figure_title = strcat(mainpath(19:end),' R-map(rlz,OQID) ',num2str(nrlz),',',OQ_RUN_ID,' Fragility ', fragility,'-',taxonomy);
risk_in = readtable(fullfile(openquakepath,strcat('damages-structural-rlz-',num2str(nrlz),'_',OQ_RUN_ID,'.csv')),'HeaderLines',1);

fault_prop = readtable(fullfile(sherifs_path1,'Faults_properties.txt'));
fault_geom = readtable(fullfile(sherifs_path1,'Faults_geometry.txt'));
fault_slip = readtable(fullfile(sherifs_path2,'mean_parameters_faults.txt'));
fault_nms = readtable(fullfile(sherifs_path2,'slip_rep_on_faults_all_data.txt'));

pre = 1:2:size(risk_in,1);
post = 2:2:size(risk_in,1);
x = risk_in.Var3(pre);
y = risk_in.Var4(pre);
if strcmp(taxonomy,'post2001')
    z = risk_in.Var10(post);
    fprintf('Doing post2001\n');
    else
    z = risk_in.Var10(pre);
    fprintf('Doing pre1981\n');
end

figure(1)

 hold on
 ax = usamap(latlim,lonlim);
 setm(ax,'MapProjection','mercator', 'MapLatLimit',latlim,'MapLonLimit',lonlim,'FlineWidth',0.7,'FontSize',6);
%

hs=scatterm(y,x,20,z,'filled');
 hs.Children.MarkerFaceAlpha = .5;
 hs.Children.MarkerEdgeColor = 'k';
 colormap(mycolors)
 risk_max = ceil(max(z)*10000)/10000;
 caxis([0 risk_max]);
 c_tick = linspace(0,risk_max,size(mycolors,1)+1);
 colbprop = colorbar;
 colbprop.Ticks = c_tick ;
 colbprop.Position = [0.60 0.60 0.014 0.24];
 colbprop.FontSize = 6;
 ylabel(colbprop,'APOC')

 % plot faulty sections
 for i = 1:size(fault_prop,1)
   %  slip = 1;
   slip_input = fault_slip.Var4(i);
   nms_output = fault_nms.Var22(i);
   slip_output = slip_input - slip_input*nms_output/100;
   slip_pos = find(histcounts(slip_output,limitisliprate));  
id = fault_prop.Var2{i};
pos = find(strcmp(fault_geom.Var2,id)==1);
xf = fault_geom.Var3(pos);
yf = fault_geom.Var4(pos);

plotm(yf,xf,'-k','LineWidth',1,'color',coloreslip(slip_pos,:))
title(figure_title,'FontSize',5, 'Interpreter', 'none');

 end
% legend of SLIPRATES

for srlim = 1:size(coloreslip,1)
srleg(srlim,1) = plotm (10,90,'-','LineWidth',1,'color',coloreslip(srlim,:),...
    'Display',char(labelsliprate(srlim)));
end
legend1 = legend(srleg);
legend1.Position = [0.38 0.25 0.07 0.12];
legend1.Box = 'off';
legend1.FontSize = 5;
title(legend1,'slip rate (mm/yr)')
legend1.Title.Visible = 'on';
hold off 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
saveas(1,fullfile(model_output,'figure',strcat('map_risk_',fragility,'_',taxonomy,'_',date,'.png')),'tiff')


