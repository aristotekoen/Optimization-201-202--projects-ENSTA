%%%  EQUILIBRE D'UNE CHAINE ARTICULEE  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Le programme suivant permet de simuler une chaine formee de barres rigides
%  contenue dans un plan vertical et de trouver sa position d'equilibre statique
%  en utilisant un algorithme d'optimisation pour minimiser son energie
%  potentielle sous contraintes appropriees.
%
%  Ce programme est constitue des fichiers suivants :
%    - ch.m   code principal contenant les procedures appelantes
%    - sqp.m  code d'optimisation : methode de Newton avec recherche lineaire
%    - chs.m  simulateur
%
%  L'utilisateur pourra adapter le programme au probleme qui l'interesse en 
%  modifiant les parametres suivants : 
%    - options.tol     permet de fixer les criteres d'arret de l'algorithme de 
%                      Newton : il s'agit des normes infinies du gradient du 
%                      Lagrangien et du vecteur des contraintes, respectivement
%                      aux positions 1 et 2
%    - options.maxit   nombre maximal d'iterations de l'aglorithme
%    - options.rl      permet d'activer la recherche lineaire (=0 <=> active)
%    - options.verb    gestion de l'affichage des variables : 
%                         =0 pas d'affichage
%                         =1 affichage des variables � chaque iteration
%                         =2 detail de la recherche lineaire
%    - run_cas_test    un certain nombre de cas-tests sont disponibles : tester 
%                      '2a', '2b', '2c', '2d', '3a', '3b', '3d', '4a'
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-------------------------------------------------------------------------------
% Parametres modifiables

options.tol = [10^(-6) 10^(-6) 10^(-6)];    % criteres d'arret
options.maxit = 100;                 % nombre max d'iterations
options.rl = 1;                     % recherche lineaire (0 = active)
options.verb = 1;                   % affichage (voir ci-dessus)
options.sqp = 0;                    % Algo SQP
options.deriv = 1;                  % algo quasi-Newton (2 = Newton)
run_cas_test = '6c';                % cas-test a appeler

%-------------------------------------------------------------------------------
% Procedures d'appel

% variables globales
global L;     % longueur des barres
global nb;    % nombre de barres
global nn;    % nombre de noeuds
global A;     % abscisse du second point de fixation
global B;     % ordonnee du second point de fixation
global R;     % vecteur R pour le plancher
global S;     % vecteur S du plancher  
global p;     % dimension de R et S
global FIG_YMIN;
global FIG_YMAX;
global FIG_XMIN;
global FIG_XMAX;


fprintf('optimisation sur le cas-test n�%s \n', run_cas_test)

% Cas-tests TP 2
if run_cas_test(1) == '2'
  
    L = [.7 .5 .3 .2 .5]';
    nb = length(L);  
    nn = nb-1;
    A=1;
    B=-1;

    lm = [1 1 1 1 1]'; % vecteur initial des multiplicateurs de Lagrange
    lmi ='no';
    % cas-test 2a 
    if run_cas_test == '2a'
      
        % position initiale des noeuds
        xy = [.2 .4 .6 .8 ...
            -1.0 -1.5 -1.5 -1.3]';
        [x, lm,lmi, info] = sqp(@chs, xy, lm,lmi, options);    % appel de l'optimiseur
    end

    % cas-test 2b
    if run_cas_test == '2b'
    
        xy = [.2 .4 .6 .8 ...
             1.0 1.5 1.5 1.3]';
        [x, lm,lmi, info] = sqp(@chs, xy, lm,lmi, options);
    end
    
    % cas-test 2c
    if run_cas_test == '2c'
      
        xy = [.2 .4 .6 .8 ...
             -1.0 -1.5 1.5 -1.3]';
        [x, lm,lmi, info] = sqp(@chs, xy, lm,lmi, options);
    end

    % cas-test 2d
    if run_cas_test == '2d'
    
        xy = [.2 .4 .6 .8 ...
              1.0 -1.2 1.5 -1.3]';
        [x, lm,lmi, info] = sqp(@chs, xy, lm,lmi, options);
    end
end

% Cas-test 3a : 
if run_cas_test == '3a'
  
    L = [0.6 0.6]';
    nb = length(L);
    nn = nb-1;
    A = 1;
    B = 0;

    xy = [0.5 0.4]';
    lm = [1 1]';
    lmi = 'no';
    [x, lm,lmi, info] = sqp(@chs, xy, lm,lmi, options);
end

% cas-test 3b :
if run_cas_test == '3b'
  
    L = [0.5 2 0.5]';
    nb = length(L);
    nn = nb-1;
    A = 1;
    B = 0;

    lm = [1 1 1]';
    lmi ='no';
    xy = [0.3 0.7 -0.2 -0.2]';
    [x, lm,lmi, info] = sqp(@chs, xy, lm,lmi, options);
    figure(1);
    xlim([-1 2]);
    ylim([-0.5 0.5]);
end

% cas-test 3c :
if run_cas_test == '3c'
  
    L = [0.5 2 0.5]';
    nb = length(L);
    nn = nb-1;
    A = 0;
    B = -1;

    lm = [1 1 1]';
    lmi ='no';
    xy = [0.5 0.5 0 -1]';
    [x, lm,lmi, info] = sqp(@chs, xy, lm,lmi, options);

    figure(1);
    xlim([-0.5 0.5]);
    ylim([-2 1]);
end

% cas-test 4a :
if run_cas_test == '4a'
  
    L = [0.7 0.5 0.3 0.2 0.5]';
    nb = length(L);
    nn = nb-1;
    A = 1;
    B = -1;

    xy = [0.2 0.4 0.6 0.8 ...
          1 1.5 1.5 1.3 ]';
    lme = [1 1 1 1 1]';
    lmi = 'no';
    p=0;
    R=0;
    S=0;
    
    [x, lm, info] = sqp(@chs, xy, lme,lmi, options);
end

if run_cas_test == '4b'
  
    L = [0.2 0.2 0.2 0.3 0.3 0.5 0.2 0.2 0.3 0.1]';
    nb = length(L);
    nn = nb-1;
    A = 1;
    B = 0;

    xy = [0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 ...
          -0.5 -0.9 -1.2 -1.4 -1.5 -1.4 -1.2 -0.9 -0.5 ]';
    lm = [0 0 0 0 0 0 0 0 0 0]';
    lmi = [0 0 0 0 0 0 0 0 0]';
    R = -0.25;
    S = -0.5;
    p = length(R);
    
    [x, lme, lmi, info] = sqp(@chs, xy, lm,lmi, options);
    figure(1);
      ylim([-0.8 0.1]);
      % dessin du plancher
      FIG_XMIN = 0;
      FIG_XMAX = 1.5;
      FIG_YMAX = 0.1;
      FIG_YMIN = -0.7;
      floor_color = 0.81*[1 1 1];
    xf = [FIG_XMIN;FIG_XMIN;FIG_XMAX;FIG_XMAX];
      for i=1:length(R)
        yf = [R(i)+S(i)*FIG_XMIN;FIG_YMIN;FIG_YMIN;R(i)+S(i)*FIG_XMAX];
        fill(xf,yf,floor_color,'EdgeColor',floor_color);
      end
  
end

if run_cas_test == '4c'
  
    L = [0.2 0.2 0.2 0.3 0.3 0.5 0.2 0.2 0.3 0.1]';
    nb = length(L);
    nn = nb-1;
    A = 1;
    B = 0;

    xy = [0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 ...
          -0.5 -0.9 -1.2 -1.4 -1.5 -1.4 -1.2 -0.9 -0.5 ]';
    lm = [1 1 1 1 1 1 1 1 1 1]';
    lmi = zeros(18,1);
    R = [-0.25 -0.5]';
    S = [-0.5 0]';
    p = length(R);
    
    
    


     
      
    [x, lme, lmii, info] = sqp(@chs, xy, lm,lmi, options);
    
    figure(1);
    ylim([-0.5 0.1]);
    % dessin du plancher
    FIG_XMIN = 0;
    FIG_XMAX = 1.5;
    FIG_YMAX = 0.1;
    FIG_YMIN = -0.5;
    floor_color = 0.81*[1 1 1];
    xf = [FIG_XMIN;FIG_XMIN;FIG_XMAX;FIG_XMAX];
    for i=1:length(R)
        yf = [R(i)+S(i)*FIG_XMIN;FIG_YMIN;FIG_YMIN;R(i)+S(i)*FIG_XMAX];
        fill(xf,yf,floor_color,'EdgeColor',floor_color);
    end
end


if run_cas_test == '5a'
   L = [.5 .75 .75 .6]';
   nb = length(L);
   nn = nb-1;
   A =1;
   B=0;
   
   
   xy = [0 .5 1 ...
         1 .75 .5]';
   R = [-1 -.5 -5]';
   S = [-4 0 4]';
   p =length(R);
   lm = [1 1 1 1]';
   lmi = zeros(9,1);
   [x, lme, lmii, info] = sqp(@chs, xy, lm,lmi, options);
   

       figure(1);
    ylim([-8 4]);
    % dessin du plancher
    FIG_XMIN = -1;
    FIG_XMAX = 2;
    FIG_YMAX = 4;
    FIG_YMIN = -8;
    floor_color = 0.81*[1 1 1];
    xf = [FIG_XMIN;FIG_XMIN;FIG_XMAX;FIG_XMAX];
    for i=1:length(R)
        yf = [R(i)+S(i)*FIG_XMIN;FIG_YMIN;FIG_YMIN;R(i)+S(i)*FIG_XMAX];
        fill(xf,yf,floor_color,'EdgeColor',floor_color);
    end
end


if run_cas_test == '6a'

   L = [3 2.5 2.5]';
   xy = [-2  0 ...
          1 -2]';
   
   R = [-6 -10]';
   S = [-2 100]';
   p = length(R);
   nb = length(L);
   nn = nb-1;
   A =0;
   B=-4;
   lm = [1 1 1 ]';
   lmi =  zeros(4,1);
   [x, lme, lmii, info] = sqp(@chs, xy, lm,lmi, options);
    
    figure(1);
    ylim([-5 1]);
    xlim([-6,2]);
    % dessin du plancher
    FIG_XMIN = -6;
    FIG_XMAX = 2;
    FIG_YMAX = 1;
    FIG_YMIN = -5;
    floor_color = 0.81*[1 1 1];
    xf = [FIG_XMIN;FIG_XMIN;FIG_XMAX;FIG_XMAX];
    for i=1:length(R)
        yf = [R(i)+S(i)*FIG_XMIN;FIG_YMIN;FIG_YMIN;R(i)+S(i)*FIG_XMAX];
        fill(xf,yf,floor_color,'EdgeColor',floor_color);
    end
end
 

if run_cas_test == '6b'

   L = [0.1 0.2 0.3 0.4 0.5 0.4 0.3 0.1]';
 %xy = [0 0.2 0.3 0.4 0.5 0.6 0.7 ...
  %       0 0 0 0 0 0 0]';
   R = [-1 -0.2 -1]';
   S = [-7 0 7]';
   p = length(R);
   nb = length(L);
   nn = nb-1;
   A =0;
   B=0;
   lm = [1 1 1 1 1 1 1 1 ]';
   lmi = zeros(21,1);
   xy = [ -0.125   -0.10   -0.25   -0.25    0.005    0.08   -0.005 ...
            0.35    0.7    1.35    2    1.3   0.6   0.3]';
   % xy = [ -0.1143   -0.1052   -0.2537   -0.2520    0.0050    0.0779   -0.0044 ...
    %        0.3496    0.6681    1.3492    1.9702    1.2846    0.6096    0.2781]';
    
    [x, lme, lmii, info] = sqp(@chs, xy, lm,lmi, options);
    
    figure(1);
    ylim([-2 2]);
    % dessin du plancher
    FIG_XMIN = -2;
    FIG_XMAX = 2;
    FIG_YMAX = 2;
    FIG_YMIN = -2;
    floor_color = 0.81*[1 1 1];
    xf = [FIG_XMIN;FIG_XMIN;FIG_XMAX;FIG_XMAX];
    for i=1:length(R)
        yf = [R(i)+S(i)*FIG_XMIN;FIG_YMIN;FIG_YMIN;R(i)+S(i)*FIG_XMAX];
        fill(xf,yf,floor_color,'EdgeColor',floor_color);
    end   
    
end

if run_cas_test == '6c'
 zoom = sqrt(0.01 + 0.01);
   L = zoom*ones(12,1);
 xy = [ 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1 1.1 ...
         -0.1 0.1 -0.1 0.1 -0.1 0.1 -0.1 0.1 -0.1 0.1 -0.1]';
   R = [-1 -0]';
   S = [2 -.5 ]'
   p = length(R);
   nb = length(L);
   nn = nb-1;
   A =0.7;
   B=0.7;
   lm = [1 1 1 1 1 1 1 1 1 1 1 1]';
   lmi = zeros(22,1);
   %xy = [ -0.1143   -0.1052   -0.2537   -0.2520    0.0050    0.0779   -0.0044 ...
    %        0.3496    0.6681    1.3492    1.9702    1.2846    0.6096    0.2781]';
   % xy = [ -0.1143   -0.1052   -0.2537   -0.2520    0.0050    0.0779   -0.0044 ...
    %        0.3496    0.6681    1.3492    1.9702    1.2846    0.6096    0.2781]';
    
    [x, lme, lmii, info] = sqp(@chs, xy, lm,lmi, options);
    
    figure(1);
    ylim([-2 2]);
    % dessin du plancher
    FIG_XMIN = -2;
    FIG_XMAX = 2;
    FIG_YMAX = 2;
    FIG_YMIN = -2;
    floor_color = 0.81*[1 1 1];
    xf = [FIG_XMIN;FIG_XMIN;FIG_XMAX;FIG_XMAX];
    for i=1:length(R)
        yf = [R(i)+S(i)*FIG_XMIN;FIG_YMIN;FIG_YMIN;R(i)+S(i)*FIG_XMAX];
        fill(xf,yf,floor_color,'EdgeColor',floor_color);
    end   
    
end