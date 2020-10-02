    clear all
    clc

   RC1 =  [-0.17	2.80; 0.64	1.42;  0.38	1.08; 0.23	0.85];
   RC2 =  [-1.34 2.06; -0.18 0.98; -0.26 0.64; 0.31 0.91; 0.67 0.93];
   RMA2 = [-10.76 9.84; 1.30 9.89; 5.48 9.33; 5.42 6.18; 3.53 3.41];
   IMA7 = [-3.64	2.89;-1.48	2.59; -0.75	1.86; 0.67	2.16];
   IMA2 = [-9.96 5.95; -3.02 5.30; -1.12 4.77; -0.21 2.95; 0.04 1.84];
   IMA1 = [-10.48 8.01; 0.04 10.81; 2.13 4.99; 3.06 3.81];
 
   fragpre1919 = [-3.10 2.89; -1.29 2.89; -0.52 2.89; 0.64 2.89; 2.43 2.89];
   frag19821991 = [-0.12 2.99;1.98 2.99;2.72 2.99; 3.61 2.99;5.00 2.99];
   fragpost2001 = [ 0.39 3.24; 1.71 3.24; 2.48 3.24; 3.45 3.24; 5.03 3.24];
 
    pre1919rota = [0.154,0.856;0.228,0.856;0.276,0.856; 0.377, 0.856;0.656, 0.856];
    pre1919rota = [log(pre1919rota(:,1)),pre1919rota(:,2)];
   
    post1981rota = [0.626,1.175;1.292,1.175;1.560,1.175; 2.081, 1.175;3.809, 1.175];
    post1981rota = [log(post1981rota(:,1)),post1981rota(:,2)];
   
   nodamage = 0.05;
   statelim={'DS1','DS2','DS3','DS4','DS5'};
 IML = load('IML.txt');
    modelli = {'pre1919rota','post1981rota'}
    %modelli = {'pre1919','post2001'}
   parametri = cat(3,pre1919rota,post1981rota);

 fidOQ=fopen('structural_fragility_model_rota.xml','w');
    fprintf(fidOQ, '<?xml version="1.0" encoding="UTF-8"?>\n\n');
    fprintf(fidOQ,'<nrml xmlns="http://openquake.org/xmlns/nrml/0.5">\n');
    fprintf(fidOQ,'<fragilityModel id="fragility_example"\n');
    fprintf(fidOQ,'assetCategory="buildings"\n');
    fprintf(fidOQ,'lossCategory="structural">\n');
    
    fprintf(fidOQ,' <description>Fragility Model Example</description>\n');
    fprintf(fidOQ,'<limitStates>DS1 DS2 DS3 DS4 DS5</limitStates>\n');
      for f = 1:size(parametri,3)
 
       fprintf(fidOQ,strcat('<fragilityFunction id="',char(modelli(f)),'" format="discrete">\n'));
  fprintf(fidOQ,strcat('<imls imt="PGA" noDamageLimit="0.05">',num2str(IML),'</imls>\n'));

         for d = 1:size(parametri,1) % 
       frag(d,:,f)  = logncdf(IML,parametri(d,1,f),parametri(d,2,f));
       fprintf(fidOQ,strcat('<poes ls="',char(statelim(d)),'">',num2str(frag(d,:,f)),'</poes>\n'));
         end
       
         fprintf(fidOQ,'</fragilityFunction>\n');
   end
fprintf(fidOQ,'</fragilityModel>\n');

fprintf(fidOQ,'</nrml>\n');

