%xml2text
clear all
clc


filename = 'source_model_A4_m1_maxmag1_storica' ;
xml=xml2struct(strcat(filename,'.xml'));
numero_sorgenti=size(xml.nrml.sourceModel.areaSource,2);


% pre allocation of matrix, note that lengths of not-NaN occrrences and
% vertexes that will be allocated are different

occurrences(1:numero_sorgenti,1:45)=NaN;
geometry(1:numero_sorgenti,1:200)=NaN;



for i=1:numero_sorgenti
    id(i,1)=str2num(xml.nrml.sourceModel.areaSource{i}.Attributes.id);
    temp_occ=[];
    temp_vertex=[];
    temp_occ=str2num(xml.nrml.sourceModel.areaSource{i}.incrementalMFD.occurRates.Text);
    occurrences(i,1:length(temp_occ))=temp_occ;
    temp_vertex=str2num(xml.nrml.sourceModel.areaSource{i}.areaGeometry.Polygon.exterior.LinearRing.posList.Text);
    geometry(i,1:length(temp_vertex))=temp_vertex;
end

% output
save(strcat(filename,'_id.txt'),'id','-ascii')
save(strcat(filename,'_occurrences.txt'),'occurrences','-ascii')
save(strcat(filename,'_geometry.txt'),'geometry','-ascii')
