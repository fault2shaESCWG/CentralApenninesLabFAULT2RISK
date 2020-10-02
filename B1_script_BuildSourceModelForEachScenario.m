%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this script creates the xml for each rupture scenario
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
clc
close all
warning('off','all')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% USER OPTIONS
% path for the SHERIFS output xml source model
% we here use only the 1st MonteCarlo run based on the mean parameters
pathin = fullfile('A_SHERIFS_CAD','CAD_optionA1B1C1','ModelMultiFaultsTest','bg_BG_1','Le2010_A_a','sc_Set_1','bmin_0.8_bmax_1.1','MFD_double_GR')
filename = 'Source_model_1';           % xml file
% make output directory for xml files
mainpath = 'WORKING_DIRECTORY_A1B1C1_10km'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pathout1 = fullfile(mainpath,'OQinputs','Models');
if isdir(pathout1)==0
mkdir (pathout1)
elseif isdir(pathout1)==1
    rmdir (pathout1,'s')
    mkdir (pathout1)
end

if isdir(pathout1)==0
mkdir (pathout1)
end
% set output input directory for OQ files
pathout= fullfile(mainpath,'OQinputs')



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% COPY and LOAD INPUTS XML
copyfile (fullfile(pathin,strcat(filename,'.xml')),fullfile(pathout1,strcat(filename,'.xml')));
xml = xml2struct( fullfile(pathin,strcat(filename , '.xml'))) ;                % xml file

numero_sorgenti_faglie = size( xml.nrml.sourceModel.simpleFaultSource,2) ;  %
numero_sorgenti_ch = size( xml.nrml.sourceModel.characteristicFaultSource,2)   %



% read simple faults sources

for i = 1:numero_sorgenti_faglie
        ibranch = i-1; 
  if ibranch <10
    branchID = strcat('br000',num2str(ibranch));
      elseif ibranch >=10 & i < 100
    branchID = strcat('br00',num2str(ibranch));
      elseif ibranch >=100 & i < 1000
    branchID = strcat('br0',num2str(ibranch));
    elseif ibranch >=1000 & i < 10000
    branchID = strcat('br',num2str(ibranch));
  end  
    name = (xml.nrml.sourceModel.simpleFaultSource{i}.Attributes.name);
    tabella(i,:) = cell2table({branchID,name});
    
    modelname=(fullfile(pathout1,strcat(filename,'_',branchID,'.xml')));
    fidOQ=fopen(modelname,'w');

fprintf(fidOQ, '<?xml version="1.0" encoding="UTF-8"?>\n\n');
fprintf(fidOQ,'<nrml xmlns="http://openquake.org/xmlns/nrml/0.4" xmlns:gml="http://www.opengis.net/gml">\n');
fprintf(fidOQ,strcat(strcat([blanks(3),'<sourceModel name="',filename,'_',date,'">\n'])));
id_source=0; % sources ID are read from ID = 0 
  
    
    
    minmag = str2num(xml.nrml.sourceModel.simpleFaultSource{i}.incrementalMFD.Attributes.minMag);
    occurrences =  (xml.nrml.sourceModel.simpleFaultSource{i}.incrementalMFD.occurRates.Text);
    rake = str2num(xml.nrml.sourceModel.simpleFaultSource{i}.rake.Text);
    dip = str2num(xml.nrml.sourceModel.simpleFaultSource{i}.simpleFaultGeometry.dip.Text); 
    coord = str2num(xml.nrml.sourceModel.simpleFaultSource{i}.simpleFaultGeometry.gml_colon_LineString.gml_colon_posList.Text);
    name = (xml.nrml.sourceModel.simpleFaultSource{i}.Attributes.name);
    id = (xml.nrml.sourceModel.simpleFaultSource{i}.Attributes.id);
    usd =     str2num(xml.nrml.sourceModel.simpleFaultSource{i}.simpleFaultGeometry.upperSeismoDepth.Text);
    lsd =     str2num(xml.nrml.sourceModel.simpleFaultSource{i}.simpleFaultGeometry.lowerSeismoDepth.Text);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
id_source=id_source+1;
fprintf(fidOQ,strcat([blanks(0),'<simpleFaultSource id="',num2str(id),'" name="',name,'" tectonicRegion="Active_Shallow_Crust">\n']));
fprintf(fidOQ,strcat([blanks(5),'<simpleFaultGeometry>\n']));
fprintf(fidOQ,strcat([blanks(7),'<gml:LineString>\n']));
fprintf(fidOQ,strcat([blanks(7),'<gml:posList>\n']));
for i_coord = 1:size(coord,1)
fprintf(fidOQ,strcat([blanks(7),num2str(coord(i_coord,:)),'\n']));
end
fprintf(fidOQ,strcat([blanks(7),'</gml:posList>\n']));
fprintf(fidOQ,strcat([blanks(7),'</gml:LineString>\n']));
fprintf(fidOQ,strcat([blanks(7),'<dip>',num2str(dip),'</dip>\n']));
    fprintf(fidOQ,strcat([blanks(7),'<upperSeismoDepth>',num2str(usd),'</upperSeismoDepth>\n']));
    fprintf(fidOQ,strcat([blanks(7),'<lowerSeismoDepth>',num2str(lsd),'</lowerSeismoDepth>\n']));
    fprintf(fidOQ,strcat([blanks(5),'</simpleFaultGeometry>\n']));
fprintf(fidOQ,strcat([blanks(7),'<rake>',num2str(rake),'</rake>\n']));
fprintf(fidOQ,strcat([blanks(7),'<magScaleRel>WC1994</magScaleRel>\n']));
fprintf(fidOQ,strcat([blanks(7),'<ruptAspectRatio>1.0</ruptAspectRatio>\n']));

fprintf(fidOQ,strcat([blanks(7),'<incrementalMFD minMag="',num2str(minmag),'" binWidth="',num2str(0.1),'">\n']));
fprintf(fidOQ,strcat([blanks(9),'<occurRates>',occurrences,'</occurRates>\n']));
fprintf(fidOQ,strcat([blanks(7),'</incrementalMFD>\n']));

fprintf(fidOQ,strcat([blanks(0),'</simpleFaultSource>\n']));


fprintf(fidOQ,strcat([blanks(3),'</sourceModel>\n']));
    fprintf(fidOQ,strcat([blanks(3),'</nrml>']));
fclose(fidOQ);
end


% reading multi-fault sources
previous_i = i-1;

for i = 1:numero_sorgenti_ch
   if (i+previous_i) <10
    branchID = strcat('br000',num2str(i+previous_i));
      elseif (i+previous_i) >=10 & (i+previous_i) < 100
    branchID = strcat('br00',num2str(i+previous_i));
      elseif (i+previous_i) >=100 & (i+previous_i) < 1000
    branchID = strcat('br0',num2str(i+previous_i));
    elseif (i+previous_i) >=1000 & (i+previous_i) < 10000
    branchID = strcat('br',num2str(i+previous_i));
   end   
      
   name = (xml.nrml.sourceModel.characteristicFaultSource{i}.Attributes.name);
tabella((i+previous_i+1),:) = cell2table({branchID,name});

  modelname=(strcat(pathout1,'/Source_model_1_',branchID,'.xml'));
fidOQ=fopen(modelname,'w');

fprintf(fidOQ, '<?xml version="1.0" encoding="UTF-8"?>\n\n');
fprintf(fidOQ,'<nrml xmlns="http://openquake.org/xmlns/nrml/0.4" xmlns:gml="http://www.opengis.net/gml">\n');
fprintf(fidOQ,strcat(strcat([blanks(3),'<sourceModel name="20200520">\n'])));
id_source=0; % sources ID are read from ID = 0 
      
    minmag = str2num(xml.nrml.sourceModel.characteristicFaultSource{i}.incrementalMFD.Attributes.minMag);
    occurrences =  (xml.nrml.sourceModel.characteristicFaultSource{i}.incrementalMFD.occurRates.Text);
    name = (xml.nrml.sourceModel.characteristicFaultSource{i}.Attributes.name);
    id = (xml.nrml.sourceModel.characteristicFaultSource{i}.Attributes.id);
    n_fault = size(xml.nrml.sourceModel.characteristicFaultSource{i}.surface.simpleFaultGeometry,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
id_source=id_source+1;
fprintf(fidOQ,strcat([blanks(0),'<characteristicFaultSource id="',num2str(id),'" name="',name,'" tectonicRegion="Active_Shallow_Crust">\n']));
fprintf(fidOQ,strcat([blanks(5),'<surface>\n']));
for nf = 1:n_fault
     coord = str2num(xml.nrml.sourceModel.characteristicFaultSource{i}.surface.simpleFaultGeometry{nf}.gml_colon_LineString.gml_colon_posList.Text);
     dip = str2num(xml.nrml.sourceModel.characteristicFaultSource{i}.surface.simpleFaultGeometry{nf}.dip.Text);
     usd = str2num(xml.nrml.sourceModel.characteristicFaultSource{i}.surface.simpleFaultGeometry{nf}.upperSeismoDepth.Text);
     lsd = str2num(xml.nrml.sourceModel.characteristicFaultSource{i}.surface.simpleFaultGeometry{nf}.lowerSeismoDepth.Text);

 fprintf(fidOQ,strcat([blanks(5),'<simpleFaultGeometry>\n']));
fprintf(fidOQ,strcat([blanks(7),'<gml:LineString>\n']));
fprintf(fidOQ,strcat([blanks(7),'<gml:posList>\n']));
for i_coord = 1:size(coord,1)
fprintf(fidOQ,strcat([blanks(7),num2str(coord(i_coord,:)),'\n']));
end
fprintf(fidOQ,strcat([blanks(7),'</gml:posList>\n']));
fprintf(fidOQ,strcat([blanks(7),'</gml:LineString>\n']));
fprintf(fidOQ,strcat([blanks(7),'<dip>',num2str(dip),'</dip>\n']));
    fprintf(fidOQ,strcat([blanks(7),'<upperSeismoDepth>',num2str(usd),'</upperSeismoDepth>\n']));
    fprintf(fidOQ,strcat([blanks(7),'<lowerSeismoDepth>',num2str(lsd),'</lowerSeismoDepth>\n']));
    fprintf(fidOQ,strcat([blanks(5),'</simpleFaultGeometry>\n']));

end
fprintf(fidOQ,strcat([blanks(5),'</surface>\n']));
fprintf(fidOQ,strcat([blanks(7),'<rake>',num2str(rake),'</rake>\n']));
fprintf(fidOQ,strcat([blanks(7),'<magScaleRel>WC1994</magScaleRel>\n']));
fprintf(fidOQ,strcat([blanks(7),'<ruptAspectRatio>1.0</ruptAspectRatio>\n']));

fprintf(fidOQ,strcat([blanks(7),'<incrementalMFD minMag="',num2str(minmag),'" binWidth="',num2str(0.1),'">\n']));
fprintf(fidOQ,strcat([blanks(9),'<occurRates>',occurrences,'</occurRates>\n']));
fprintf(fidOQ,strcat([blanks(7),'</incrementalMFD>\n']));

fprintf(fidOQ,strcat([blanks(0),'</characteristicFaultSource>\n']));


fprintf(fidOQ,strcat([blanks(3),'</sourceModel>\n']));
    fprintf(fidOQ,strcat([blanks(3),'</nrml>']));
fclose(fidOQ);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% building SourceModelLogicTree

scenarios = numel(dir(fullfile(pathout1,'*.xml')));
scenarios = scenarios -1;

n_branches = scenarios +1 % add the source model "total" to the count of the number of sources

pesi(1:n_branches,1)=round(1/n_branches,7);
pesi(end,1) = 1-sum(pesi(1:(end-1),1));

fidOQ=fopen(fullfile(pathout,'Sources_Logic_tree.xml'),'w');
fprintf(fidOQ,'<?xml version="1.0" encoding="UTF-8"?>\n');
fprintf(fidOQ,'<nrml xmlns:gml="http://www.opengis.net/gml"\n');
fprintf(fidOQ,'      xmlns="http://openquake.org/xmlns/nrml/0.4">\n' );
fprintf(fidOQ,'<logicTree logicTreeID="lt1">\n' );
fprintf(fidOQ,'   <logicTreeBranchingLevel branchingLevelID="bl1">\n' );
fprintf(fidOQ,'        <logicTreeBranchSet uncertaintyType="sourceModel" branchSetID="bs1">\n' );

 for i = 1:(n_branches-1)
     branch = i-1;
      if branch <10
    branchID = strcat('br000',num2str(branch));
      elseif branch >=10 & branch < 100
    branchID = strcat('br00',num2str(branch));
      elseif branch >=100 & branch < 1000
    branchID = strcat('br0',num2str(branch));
    elseif branch >=1000 & branch < 10000
    branchID = strcat('br',num2str(branch));
      end
      
 source   = strcat('Models/',filename,'_',branchID,'.xml')
    
fprintf(fidOQ,strcat([blanks(5),'<logicTreeBranch branchID="',branchID,'">\n']));
fprintf(fidOQ,strcat([blanks(7),'<uncertaintyModel>',char(source),'</uncertaintyModel>\n']));
fprintf(fidOQ,strcat([blanks(7),'<uncertaintyWeight>',num2str(pesi(i),'%.10f\t'),'</uncertaintyWeight>\n']));
fprintf(fidOQ,strcat([blanks(5),'</logicTreeBranch>\n']));
 end
 
 % add the "total" source model
 i = i+1;branch = i-1;
 branchID = strcat('br',num2str(branch));
 source   = strcat('Models/',filename,'.xml');
    
fprintf(fidOQ,strcat([blanks(5),'<logicTreeBranch branchID="',branchID,'">\n']));
fprintf(fidOQ,strcat([blanks(7),'<uncertaintyModel>',char(source),'</uncertaintyModel>\n']));
fprintf(fidOQ,strcat([blanks(7),'<uncertaintyWeight>',num2str(pesi(i),'%.10f\t'),'</uncertaintyWeight>\n']));
fprintf(fidOQ,strcat([blanks(5),'</logicTreeBranch>\n']));
%%%%


fprintf(fidOQ,'        </logicTreeBranchSet>\n');
fprintf(fidOQ,'   </logicTreeBranchingLevel>\n');
fprintf(fidOQ,' </logicTree>\n');
fprintf(fidOQ,'</nrml>\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
writetable(tabella,fullfile(mainpath,'Visualization','Data4Maps','tabella_corrispondenze.txt'))

