function [e, c,ci, g, a, ai, hl, indic] = chs(indic, xy, lm, lmi)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulateur chs.m
% 
% INPUT  - indic  pilote le comportement du simulateur :
%                   = 2  chs trace la chaine
%                   = 3  chs calcule e et c
%                   = 4  chs calcule e, c, g et a
%                   = 5  chs calcule hl
%                   = 6  chs calcule e, c, g, a, hl et trace la chaine
%        - xy     vecteur contenant la position des noeuds             (2nn x 1)
%                 (abscisses et ordonnees concatenee)
%        - lm     vecteur contenant les multiplicateurs de Lagrange     (nb x 1)
%      
% OUTPUT  - e      valeur de l'energie potentielle en xy
%         - c      vecteur des contraintes en xy                        (nb x 1)
%         - g      gradient de e en xy                                 (2nn x 1)
%         - a      matrice jacobienne des contraintes en xy           (nb x 2nn)
%         - hl     matrice hessienne du Lagrangien                   (2nn x 2nn)
%         - indic  decrit le resultat de la simulation
%                    = 0  sortie normale
%                    = 1  parametres d'entree non corrects
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %specification des variables globales
    global L;
    global nb;
    global nn;
    global A;
    global B;
    global p;
    global R;
    global S;
    
    %copie de indic dans une variable locale
    i = indic;
    
    %verification de la conformite des input
    if 2*(length(L)-1) == length(xy)
        d = sqrt(A*A+B*B);
        s = sum(L);
        if s >= d
            indic = 0;
        else
            indic = 1;
        end
    else 
        indic = 1;
    end
    
    % initialisation des variables a renvoyer
    % (si certaines variables ne sont pas demandees alors le simulateur 
    % renvera des vecteurs de taille appropriee contenant des 0)
    e = 0;
    c = zeros(nb,1);
    ci = zeros(nn*p,1);
    g = zeros(2*nn,1);
    a = zeros(nb,2*nn);
    hl = zeros(2*nn);
    ai = zeros(nn*p, 2*nn);
    
    % affichage de la chaine sur un graphique
    if i==2 || i==6
        figure(1);
        clf;
        plot([0 xy(1:nn)' A]', [0 xy(nn+1:2*nn)' B]');
        hold on
        plot([0 xy(1:nn)' A]', [0 xy(nn+1:2*nn)' B]','o');
    end
    
    % calcul de e et c
    if i==3 || i==4 || i==6
        e = 1./2 * L' * ([xy(nn+1:2*nn) ; B] + [0 ; xy(nn+1:2*nn)]);
        c = diag([xy(1:nn)' A]' - [0 xy(1:nn)']') * ([xy(1:nn)' A]' - [0 xy(1:nn)']')...
        + diag([xy(nn+1:2*nn)' B]' - [0 xy(nn+1:2*nn)']') * ([xy(nn+1:2*nn)' B]' - [0 xy(nn+1:2*nn)']')...
        - (diag(L)*L);
        if strcmp(lmi,'no')==1
            ci = [];
        else
            ci = kron(xy(1:nn),S) + kron(ones(nn,1),R) -kron(xy((nn+1):2*nn), ones(p,1));
        end
        % calcul de g et a
        if i==4 || i==6
            g = [zeros(nn,1)' 1./2*(L(1:nn) + L(2:nb))']';
            a = [spdiags([2*(xy(1:nn) - [0 xy(1:nn-1)']') ...
                         -2*([xy(2:nn)' A]' - xy(1:nn)  )],...
                         [0 -1], nb, nn) ...
                 spdiags([2*(xy(nn+1:2*nn) - [0 xy(nn+1:2*nn-1)']') ...
                         -2*([xy(nn+2:2*nn)' B]' - xy(nn+1:2*nn)  )],...
                         [0 -1], nb, nn)];
               if strcmp(lmi,'no')==1
                   ai = zeros(nn*p,2*nn);
               else
                   ai = [kron(eye(nn),S) kron(eye(nn),-ones(p,1))];
               end
            
        end
    end
   
    % calcul de hl
    if i==5 || i==6
        hl = 2*(diag([lm(1:nn)' lm(1:nn)'] + [lm(2:nb)' lm(2:nb)']) - ...
                diag([lm(2:nn)' 0 lm(2:nn)'], -1) - diag([lm(2:nn)' 0 lm(2:nn)'], 1));
    end
    
end