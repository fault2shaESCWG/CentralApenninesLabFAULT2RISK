% This script calculates the contribution of each fault section to the
% total hazard for a given afoe and for a given site.

clear all
clc
close all
warning('off','all')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% USER OPTIONS
afoe = 0.0021;
OQ_RUN_ID = '6';% Number of Openquake run ID
fprintf(['Warning: You are using OQ_RUN_ID ',num2str(OQ_RUN_ID)]);

site = [13.4 42.35]; sito = 'AQ'; % coordinate of the site
%site = [13.44 42.04]; sito = 'AV'; % coordinate of the site

mainpath = 'WORKING_DIRECTORY_A1B1C1_10km';
openquakepath = fullfile(mainpath,'OQoutputs');
model_output = fullfile(mainpath,'Visualization');
sherifs_path1 = fullfile('A_SHERIFS_CAD','data','CAD_optionA1B1C1_10km');
sherifs_path2 = fullfile('A_SHERIFS_CAD','CAD_optionA1B1C1_10km','analysis','txt_files');

limitisliprate = [0 0.1 0.5 1 3];
coloreslip = [.5 .5 .5; 0 1 0; 0 0 1; 1 0 0];
labelsliprate = {'<0.1','0.1-0.5','0.6-1.0','>1.0'};
latlim=([41.6 43.2]);
lonlim=([12.7 14.3]);

% bar plot of the contribution for INTESITY measure levels provided in the
% OQ input file : for example: 3 corresponds to the third IML value given in the
% job.ini  
levels = [3,7,17];

%threshold for identification of sections in map
thershold_participation_label_map = 95;
thershold_participation_label_figure = 3;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Openquake files: rlz corrisponding to the last file (the total hazard)
realization = readtable(fullfile(openquakepath,strcat('realizations_',num2str(OQ_RUN_ID),'.csv')));
nrlz= size(realization,1)-1; %  note that the branchID goes from 0 to n-1
HCfaulttab = readtable(fullfile(openquakepath,strcat('hazard_curve-rlz-',num2str(nrlz),'-PGA_',OQ_RUN_ID,'.csv')),'HeaderLines',1);
iml = HCfaulttab.Properties.VariableNames;
IML_1=strrep(iml,'poe_','');IML_1=strrep(IML_1,'_','.');
IML_1 = (IML_1(4:end));
IML=str2num(char(IML_1));
IML = IML';

figure_title = strcat(model_output(19:29),'-',sito,'-',num2str(nrlz),'-',OQ_RUN_ID);
%%%%%%%%%%%%
% inputs
modelli = readtable(fullfile(model_output,'Data4Maps','tabella_corrispondenze.txt'),'Delimiter',',');
fault_prop = readtable(fullfile(sherifs_path1,'Faults_properties.txt'));
fault_geom = readtable(fullfile(sherifs_path1,'Faults_geometry.txt'));
fault_slip = readtable(fullfile(sherifs_path2,'mean_parameters_faults.txt'));
fault_nms = readtable(fullfile(sherifs_path2,'slip_rep_on_faults_all_data.txt'));
slip_rateDP = readtable(fullfile(model_output,'Data4Maps','USED_sliprateDP.txt'));
section_name = unique(fault_geom.Var2);
IDsection = 1:size(section_name);

 % lon lat depth poe
 HCfault =table2array(HCfaulttab);

 x = HCfault(:,1);
 y = HCfault(:,2);
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% find the row in the hazard curves table corresponding to the site
 dist = sqrt( (HCfaulttab.lon-site(1)).^2 + (HCfaulttab.lat-site(2)).^2);
 s = find(dist == min(dist),1,'first');
 HC = table2array(HCfaulttab(s,4:end));
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% prepare the table with the HC at the site for the n-scenario
for j =1:nrlz % read the HC for each scenario
    i=j-1;
    T = [];
    if i>=100
        T = readtable(fullfile(openquakepath,strcat('hazard_curve-rlz-',num2str(i),'-PGA_',OQ_RUN_ID,'.csv')),'HeaderLines',1);
    elseif i <100 & i>=10
        T = readtable(fullfile(openquakepath,strcat('hazard_curve-rlz-0',num2str(i),'-PGA_',OQ_RUN_ID,'.csv')),'HeaderLines',1);    
    elseif i <10
         T = readtable(fullfile(openquakepath,strcat('hazard_curve-rlz-00',num2str(i),'-PGA_',OQ_RUN_ID,'.csv')),'HeaderLines',1);
    end
     
   tabella(j,:) = [i, table2array(T(s,4:end))];
    
end  

% for each section find the scenario that contains it and sum the HC
% then normalize the total HC for each section to ensure that the sum of HC
% of each scenario is the total HC
% this is called HCSections
% contributions in percentage are given by PartecipationSections

SommaHazard=[];pos=[];
for ns = 1:size(section_name,1)
section = (section_name{ns}) 
[u]=strfind(modelli.Var2,section);
for i = 1:nrlz
    pos(i,ns) = ~isempty(u{i});
end
temp = [];
temp = sum(tabella(find(pos(:,ns)>0),2:end),1);
SommaHazard(ns,:) = (temp);
end
temp2 = sum(SommaHazard,1);
temp3 = SommaHazard./repmat(temp2,ns,1);

HCSections = temp3 .*repmat(HC,ns,1);
PartecipationSections = HCSections./repmat(HC,ns,1);
SumHC = sum(HCSections,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot the scenario hazard curves
figure(1)
hold on
plot(IML,HC,'-r') % the total curve at the site
plot(IML,HCSections) % sections HCs
plot(IML,SumHC,'--k') % sum of HCsections should be equal to HC
title(figure_title, 'Interpreter', 'none');
hold off
set (gca, 'YScale','log','XScale','log')
xlim([IML(1) IML(end)])
grid on
saveas(1,fullfile(model_output,'figure',strcat(char(sito),'_HC_',date,'.png')),'tiff')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


figure(2)
hold on
for l = 1: length(levels)
subplot(1,length(levels),l)
hold on
bar(PartecipationSections(:,levels(l))*100);
for i = 1:size(PartecipationSections,1)
    if PartecipationSections(i,levels(l)) > 0.02
        testo=(section_name{i,:});
        testo = strrep(testo,'_','');
        text(i,PartecipationSections(i,levels(l))*100,testo,'Rotation',90,'FontSize',4)
    end
end
hold off
xlabel('sections');
ylabel(['percentage of contribution to IML:', num2str(IML(levels(l)))])
title(figure_title, 'Interpreter', 'none');
end
saveas(2,fullfile(model_output,'figure',strcat(char(sito),'_barplot_',date,'.png')),'tiff')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%imagesc plot of the contribution for some levels 
%% 
figure(3)
hold on
% Create axes
% axes1 = axes('Parent',3,...
%     'Position',[0.1 0.11 0.754792626728111 0.8]);
% hold(axes1,'on');

imagesc(PartecipationSections*100);
%xtickformat('%3.2f');
set(gca,'XGrid','off','YGrid','on','GridLineStyle','-',...
    'XTick',2:2:19,'XTickLabel',{round(IML(2:2:19),2)},'XTickLabelRotation',90,...
    'YTick',5:5:82,'YTickLabel',{5:5:82},'FontSize',10)

xlim([0.5 48.5])
ylim([0 83.5])
colormap(cool)
cbar = (colorbar)
%cbar.Location = ('southoutside');
%cbar.Position = ([0.120737327536606 0.890818858560795 0.774193548039431 0.0235732152288829]);
%cbar.Position =([0.1298    0.2095    0.7750    0.0208]);
ylabel(cbar,'Section Participation to Total Hazard (%)','FontSize', 12)
xlabel('IML (PGA)','FontSize', 12);
ylabel('section','FontSize', 12)


line([7 7],[0 83.5])
line([17 17],[0 83.5])
% add text with id section at the end
xpos1 = repmat(20,90,1);
xpos2 = repmat([31;21],45,1);
ipos = 1
for i =1:size(PartecipationSections,1)
if sum((PartecipationSections(i,:)) > thershold_participation_label_figure/100)>0
        testo=(section_name{i,:});
        testo = strrep(testo,'_','');
        text(xpos2(ipos),i,testo,'FontSize',6)
        line([xpos1(i) xpos2(ipos)],[i i])
        ipos = ipos+1;
end
end

%title(figure_title, 'Interpreter', 'none');
saveas(3,fullfile(model_output,'figure',strcat(char(sito),'_imagesc_',date,'.png')),'tiff')
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%

for l=1:length(levels)
yl = 43.00; % to place the legend
xl = 13.33; % to place the legend

figure(4)
hold on
ax = usamap(latlim,lonlim);
setm(ax,'MapProjection','mercator', 'MapLatLimit',latlim,'MapLonLimit',lonlim,'FlineWidth',0.7,'FontSize',6);
textm(yl,xl,strcat('IML "PGA" : ',num2str(IML(levels(l)),'%4.3f'),'(g)'),'FontSize',6,'FontAngle','Italic');
for ns = 1:size(PartecipationSections,1)
section = (section_name{ns}) 
pos1 = find((strcmpi(fault_slip.Var3,section))==1);
pos2 = find((strcmpi(fault_geom.Var2,section))==1);

slip_input = fault_slip.Var4(pos1);
nms_output = fault_nms.Var22(pos1);
slip_output = slip_input - slip_input*nms_output/100;
slip_pos = find(histcounts(slip_output,limitisliprate)); 
section_coord = [fault_geom.Var3(pos2),fault_geom.Var4(pos2)];
section_z = repmat(0,length(pos2),1);
hs = [];
hs = plot3m(section_coord(:,2),section_coord(:,1),section_z,'-');
hs.LineWidth =min([80*(PartecipationSections(ns,levels(l)))+0.01,3]);
hs.Color = coloreslip(slip_pos,:);

pct_level= prctile(PartecipationSections(:,levels(l)),thershold_participation_label_map);

if PartecipationSections(ns,levels(l)) >=  pct_level
   yc = mean(section_coord(:,2))+0.02;
   xc =mean(section_coord(:,1))+0.02;
   testo1 = strcat('(',num2str(ns),')');
   textm(yc,xc,testo1,'FontSize',6,'BackgroundColor','w','Margin',0.001)
   
   testo2 = strcat(testo1,'-',section,':',num2str(PartecipationSections(ns,levels(l))*100,'%3.2f'),'%');
   testo2 = strrep(testo2,'_','');
   yl = yl-0.040;
   textm(yl,xl,testo2,'FontSize',6)
   
end


end
h = textm(site(:,2),site(:,1)-0.23,'L Aquila','FontSize',6);
t = plotm(site(:,2),site(:,1),'s','MarkerSize',8);
t.MarkerFaceColor = 'yellow';
t.MarkerEdgeColor = 'black';
title(figure_title, 'Interpreter', 'none');

% legend of SLIPRATES

for srlim = 1:size(coloreslip,1)
srleg(srlim,1) = plotm (90,90,'-','LineWidth',1,'color',coloreslip(srlim,:),...
    'Display',char(labelsliprate(srlim)));
end
legend1 = legend(srleg);
legend1.Position = [0.38 0.25 0.07 0.12];
legend1.Box = 'off';
legend1.FontSize = 6;
title(legend1,'slip rate (mm/yr)')
legend1.Title.Visible = 'on';
hold off 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hold off
saveas(4,fullfile(model_output,'figure',strcat(char(sito),'_mappa_sliprate_HCcontribution',num2str(levels(l)),'_',date,'.png')),'tiff')
close(4)
end
%%


