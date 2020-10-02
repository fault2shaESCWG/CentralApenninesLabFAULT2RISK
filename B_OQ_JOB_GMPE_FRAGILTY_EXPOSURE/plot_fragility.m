clear all
clc
% plot fragility
   RC1 =  [-0.17	2.80; 0.64	1.42;  0.38	1.08; 0.23	0.85];
   RC2 =  [-1.34 2.06; -0.18 0.98; -0.26 0.64; 0.31 0.91; 0.67 0.93];
   RMA2 = [-10.76 9.84; 1.30 9.89; 5.48 9.33; 5.42 6.18; 3.53 3.41];
   IMA7 = [-3.64	2.89;-1.48	2.59; -0.75	1.86; 0.67	2.16];
   IMA2 = [-9.96 5.95; -3.02 5.30; -1.12 4.77; -0.21 2.95; 0.04 1.84];
   IMA1 = [-10.48 8.01; 0.04 10.81; 2.13 4.99; 3.06 3.81];
   
    fragpre1919 = [-3.10 2.89; -1.29 2.89; -0.52 2.89; 0.64 2.89; 2.43 2.89];
    frag19821991 = [-0.12 2.99;1.98 2.99;2.72 2.99; 3.61 2.99;5.00 2.99];
    fragpost2001 = [ 0.39 3.24; 1.71 3.24; 2.48 3.24; 3.45 3.24; 5.03 3.24];
    frg2B5B = [-3.86, 2.05;-2.26 2.05;-1.55 2.05;-0.54 2.05;1.30 2.05];
    % rota 2020
    pre1919rota = [0.154,0.856;0.228,0.856;0.276,0.856; 0.377, 0.856;0.656, 0.856];
    pre1919rota = [log(pre1919rota(:,1)),pre1919rota(:,2)];
    post1981rota = [0.626,1.175;1.292,1.175;1.560,1.175; 2.081, 1.175;3.809, 1.175];
    post1981rota = [log(post1981rota(:,1)),post1981rota(:,2)];
   
    
    nodamage = 0.05;
   
    statelim={'DS1','DS2','DS3','DS4','DS5'};
    IML = load('IML.txt');
    %modelli = {'fragpre1919','pre1919rota'};
    modelli = {'pre1919rota','post1981rota'};
    %modelli = {'frg2B5B'};
   parametri = cat(3,pre1919rota,post1981rota)
    %parametri = cat(3,fragpre1919,pre1919rota);

    % CDF
   for d = 1:size(parametri,1) % classi di danno
       for f = 1:size(parametri,3)
       frag_cum(d,:,f)  = logncdf(IML,parametri(d,1,f),parametri(d,2,f));
       end
   end
       % PDF
   for d = 1:size(parametri,1) % classi di danno
       for f = 1:size(parametri,3)
       frag_pdf(d,:,f)  = lognpdf(IML,parametri(d,1,f),parametri(d,2,f));
       end
   end
   
   figure(1)
   hold on
   for f = 1:size(parametri,3)
   plot(IML,frag_cum(1:5,:,f)*100,'-','color',[repmat(f/(size(parametri,3)+1),1,3)],'display',char(modelli(f)))
   end
   ylabel('Probability of exceeding')
   xlabel('PGA (g)')
   grid on
   xlim([0 1.0])
   ylim([0 100])
   set(gca,'fontsize',14)
   
   figure(2)
   hold on
   for f = 1:size(parametri,3)
   plot(IML,frag_pdf(1:5,:,f)*1,'-','color',[repmat(f/(size(parametri,3)+1),1,3)],'display',char(modelli(f)))
   end
   ylabel('PDF for damage')
   xlabel('PGA (g)')
   grid on
   xlim([0 0.5])
  % ylim([0 100])
   set(gca,'fontsize',14)
   