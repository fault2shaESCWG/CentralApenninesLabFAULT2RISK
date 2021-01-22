% % This script calculates the contribution of each fault section to the
% total risk for a given damage state and for a given site 
clear all
clc
close all
warning('off','all')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% USER OPTIONS
site = [13.4 42.35]; sito = 'L Aquila'; % coordinate of the site
%site = [13.44 42.04]; sito = 'AV'; % coordinate of the site
fragility ='_Rostietal.2020-L-type';
OQ_RUN_ID = '8';% Number of Openquake run ID
fprintf(['Warning: You are using OQ_RUN_ID ',num2str(OQ_RUN_ID)]);

mainpath = 'WORKING_DIRECTORY_A1B1C1_10km';
openquakepath = fullfile(mainpath,'OQoutputs');
model_output = fullfile(mainpath,'Visualization');
sherifs_path1 = fullfile('A_SHERIFS_CAD','data','CAD_optionA1B1C1_10km');
sherifs_path2 = fullfile('A_SHERIFS_CAD','CAD_optionA1B1C1_10km','analysis','txt_files');

limitisliprate = [0 0.1 0.5 1 3];
coloreslip = [.5 .5 .5; 0 1 0; 0 0 1; 1 0 0];

%threshold for identification of sections in figures
thershold_participation_label = 90;

% latlim=([41.6 43.2]);
% lonlim=([12.7 14.3]);
latlim=([41.9 42.8]);
lonlim=([13.2 13.8]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Openquake files: rlz corresponding to the last file (the total hazard)
realization = readtable(fullfile(openquakepath,strcat('realizations_',num2str(OQ_RUN_ID),'.csv')));
nrlz= size(realization,1)-1; % note that the branchID goes from 0 to n-1

figure_title = strcat(mainpath(19:29),'-',sito,'-',num2str(nrlz),'-',OQ_RUN_ID);

%Openquake files: rlz corresponding to the last file (the total hazard)
Tot = readtable(fullfile(openquakepath,strcat('damages-structural-rlz-',num2str(nrlz),'_',OQ_RUN_ID,'.csv')),'HeaderLines',1);%A1B1C1

% inputs
modelli = readtable(fullfile(model_output,'Data4Maps','tabella_corrispondenze.txt'),'Delimiter',',');
fault_prop = readtable(fullfile(sherifs_path1,'Faults_properties.txt'));
fault_geom = readtable(fullfile(sherifs_path1,'Faults_geometry.txt'));
fault_slip = readtable(fullfile(sherifs_path2,'mean_parameters_faults.txt'));
fault_nms = readtable(fullfile(sherifs_path2,'slip_rep_on_faults_all_data.txt'));
slip_rateDP = readtable(fullfile(model_output,'Data4Maps','USED_sliprateDP.txt'));
section_name = unique(fault_geom.Var2);
IDsection = 1:size(section_name);


%find the row in the hzard curves table corresponding to the site
 dist = sqrt( (Tot.Var3-site(1)).^2 + (Tot.Var4-site(2)).^2);
 s = find(dist == min(dist),2,'first'); % two fragility
 
for j =1:(nrlz-1)
    i=j-1;
    T = [];
    if i>=100
        T = readtable(fullfile(openquakepath,strcat('damages-structural-rlz-',num2str(i),'_',OQ_RUN_ID,'.csv ')));
    elseif i <100 & i>=10
        T = readtable(fullfile(openquakepath,strcat('damages-structural-rlz-0',num2str(i),'_',OQ_RUN_ID,'.csv ')));
        
    elseif i <10
         T = readtable(fullfile(openquakepath,strcat('damages-structural-rlz-00',num2str(i),'_',OQ_RUN_ID,'.csv ')));
    end
     
   tabella(j,:) = [i, T.DS5(s(1)),T.DS5(s(2))]; % pre1919 and post2001
    
end    
%
riskofcollapse = Tot.Var10(s(1:2))
sommaGMPE1_pre1919 = sum(tabella(:,2))
sommaGMPE1_post2001 = sum(tabella(:,3))

%%
for ns = 1:size(section_name,1)
section = (section_name{ns}); 
[u]=strfind(modelli.Var2,section);
for i = 1:(nrlz-1)
    pos(i,ns) = ~isempty(u{i});
end

SommaContributi(ns,:) = table({section},sum(tabella(pos(:,ns),2)),sum(tabella(pos(:,ns),3)));
end
PartecipationSections = SommaContributi;
PartecipationSections.Var3 = (PartecipationSections{:,3}./sum(PartecipationSections{:,3}));
PartecipationSections.Var2 = (PartecipationSections{:,2}./sum(PartecipationSections{:,2}));

for f = 1:2 % 2 fragility
figure(1)
hold on
bar(PartecipationSections{:,(f+1)}*100);
for i = 1:size(PartecipationSections,1)
    if PartecipationSections{i,(f+1)} > 0.02
        testo=PartecipationSections.Var1(i);
        testo = strrep(testo,'_','');
        text(i,PartecipationSections{i,(f+1)}*100,testo,'Rotation',90)
    end
end
hold off
xlabel('sections');
ylabel('contribution in percentage to the risk of collapse')

title(figure_title, 'Interpreter', 'none');
saveas(1,fullfile(model_output,'figure',strcat(char(sito),fragility,'_barplot_riskcontribution_fragility_',num2str(f),'_',date,'.png')),'tiff')
close(1)
end

%%
 
 for f = 1:2 % 2 fragility
     
 figure(2)

 hold on
 ax = usamap(latlim,lonlim);
 ax.FontSize = 6;
 setm(ax,'MapProjection','mercator', 'MapLatLimit',latlim,'MapLonLimit',lonlim,...
     'FlineWidth',0.7,'FontSize',6,'MLabelLocation',0.2,'PLabelLocation',0.4);
setm(ax,'MapProjection','mercator', 'MapLatLimit',latlim,'MapLonLimit',lonlim,...
     'FlineWidth',0.7,'FontSize',8,'MLabelLocation',0.3,'PLabelLocation',0.3,'LabelFormat','none');


%yl = 43.032;
yl =42.699;

xl = 13.85;
for ns = 1:size(PartecipationSections,1)
section = (PartecipationSections{ns,1}) ;
pos1 = find((strcmpi(fault_slip.Var3,section))==1);
pos2 = find((strcmpi(fault_geom.Var2,section))==1);
slip_input = fault_slip.Var4(pos1);
nms_output = fault_nms.Var22(pos1);
slip_output = slip_input - slip_input*nms_output/100;
slip_pos = find(histcounts(slip_output,limitisliprate));

section_coord = [fault_geom.Var3(pos2),fault_geom.Var4(pos2)];
hs = [];
hs = plotm(section_coord(:,2),section_coord(:,1),'-');
hs.LineWidth =min([80*(PartecipationSections{ns,(f+1)}),3]);
hs.Color = coloreslip(slip_pos,:);

title(figure_title, 'Interpreter', 'none');

%threshold for identification of sections in figures
pct_level= prctile(PartecipationSections{:,2},thershold_participation_label);

if PartecipationSections{ns,2} >= pct_level
   yc = mean(section_coord(:,2))+0.02;
   xc =mean(section_coord(:,1))+0.02;
   testo1 = strcat('(',num2str(ns),')');
   textm(yc,xc,testo1,'FontSize',6,'BackgroundColor','w','Margin',0.0001)
   
    testo2 = strcat(testo1,'-',section,':',num2str(PartecipationSections{ns,2}*100,'%3.2f'),'%');
    testo2 = strrep(testo2,'_','');
    yl = yl-0.04;
    textm(yl,xl,testo2,'FontSize',6)
   
end

end
h = textm(site(:,2),site(:,1)-0.13,sito,'FontSize',6);
l= plotm(site(:,2),site(:,1),'s','MarkerSize',8);
l.MarkerFaceColor = 'yellow';
l.MarkerEdgeColor = 'black';
plotm(slip_rateDP.lat,slip_rateDP.lon,'ok','MarkerSize',4,'MarkerFaceColor','w')
title(strcat(figure_title,'fragility-',num2str(f)), 'Interpreter', 'none');

% legend of SLIPRATES
limitisliprate = [0 0.1 0.5 1 3];
coloreslip = [.5 .5 .5; 0 1 0; 0 0 1; 1 0 0];
labelsliprate = {'<0.1','0.1-0.5','0.6-1.0','>1.0'};


for srlim = 1:size(coloreslip,1)
srleg(srlim,1) = plotm (90,90,'-','LineWidth',1,'color',coloreslip(srlim,:),...
    'Display',char(labelsliprate(srlim)));
end
legend1 = legend(srleg);
%legend1.Position = [0.38 0.25 0.07 0.12];
legend1.Position = [0.53 0.73 0.07 0.12];
legend1.Box = 'off';
legend1.FontSize = 6;
title(legend1,'slip rate (mm/yr)')
legend1.Title.Visible = 'on';

plotm(90,90,'ok','MarkerSize',4,'MarkerFaceColor','w','Display','Data points');

hold off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

saveas(2,fullfile(model_output,'figure',strcat(char(sito),fragility,'_riskcontribution_fragility_',num2str(f),'_',date,'.png')),'tiff')
print(fullfile(model_output,'figure',strcat(char(sito),fragility,'_riskcontribution_fragility_',num2str(f),'_',date,'.tiff')),'-dtiff','-r600');
close(2)
 end