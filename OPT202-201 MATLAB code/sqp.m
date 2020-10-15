function [x, lme,lmi, info] = sqp(simul, x, lme,lmi, options)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Optimiseur sqp.m
% 
% INPUT  - simul    specification du simulateur (ex : @chs)
%        - x        positions initiales des noeuds                     (2nn x 1)
%                   (abscisses et ordonnees concatenee)
%        - lme      vecteur initial des multiplicateurs de              (nb x 1)
%                   Lagrange associe aux contraintes d'egalite
%        - options  parametres d'optimisation (voir details ch.m)
%      
% OUTPUT  - x      vecteur des positions finales                       (2nn x 1)
%         - lme    multiplicateurs de Lagrange associe aux              (nb x 1) 
%                  contraintes d'egalite
%         - info   decrit l'etat final de l'algorithme :
%                     * info.status = 0  terminaison normale
%                     *             = 1  INPUT inconsistants
%                     *             = 2  depassement de options.maxit
%                     * info.niter   donne le nombre d'iterations effectuees
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    global nn;
    % constantes
    omega = 10^(-4);  % intervient dans la regle d'Armijo 
    n = length(x);
    m = length(lme);
    
    % initialisation des variables
    deltaF = [];  % evaluation de la convergence etape par etape
    deltaFk =[];
    indic = 0;    % communication entre l'optimiseur et le simulateur
    flag = 0;     % etat de la condition d'optimalite
    nsim = 0;     % nombre d'appels au simulateur
    k = 0;        % compteur d'iterations
    
    if options.deriv == 1
        mi = length(lmi);
        
        % appel au simulateur
        [f,c,ci, gf, jc,ai, hl, indic] = simul(6, x, lme,lmi);
        
        if strcmp(lmi,'no')== 1
            gl = gf + jc'*lme; 
        else 
            gl = gf + jc'*lme + ai'*lmi;
        end
        M = eye(2*nn);
        norm_gl = max(abs(gl));
        norm_ce = max(abs(c));
        norm_diff = max(abs(min(lmi,-ci)));
        
        if strcmp(lmi,'no')==1
            norm_diff = 0;
        else    
            norm_diff = max(abs(min(lmi,-ci)));
        end
        if strcmp(lmi,'no')==0
            while k<= options.maxit

              if [norm_gl,norm_ce,norm_diff]<= options.tol
                  if options.verb == 1
                        fprintf('----------------------------------------------------------------------\n')
                  elseif options.verb == 2
                        fprintf('--------------------------------------------------------------\n')
                  end
                  break
              end

              %sinon on résout (4.2)
              param = optimset('Display', 'off');
              [dk, fval, exitflag,output, lambda] = quadprog(M,gf,ai,-ci,jc,-c,[],[],[],param);
              
              if (options.verb == 1) && (k == 0)
                fprintf('\n')
                fprintf('Algorithme de type quasi Newton \n')
                fprintf('----------------------------------------------\n')    
                fprintf('%4s  %11s  %11s  %14s  %8s  %8s  %8s  %8s  %8s \n',...
                        'iter', '    |gl|   ', '   |ce|    ','   (ci,lmi)   ','  |lm|  ',  '   |x|  ', ' alpha  ',' Powell ',' cond(M)')
              end
              %calcul des nouveaux x, lambda:
              alpha = 1;
              x_new = x + alpha*dk;
              lm_new = lme + alpha*(lambda.eqlin-lme);
              lmi_new = lmi + alpha*(lambda.ineqlin-lmi);
              [f_new,c_new, ci_new, gf_new, jc_new,ai_new, hl_new, indic] = simul(6, x_new, lm_new,lmi_new);
              gl_new = gf_new + jc_new'*lm_new + ai_new'*lmi_new;
              
              
              gl_var = gf + jc'*lm_new +ai'*lmi_new;
              delta_k = alpha*dk;
              vec = M*delta_k;
              gamma_kl = gl_new -gl_var;


                  if gamma_kl'*delta_k < 0.2*(delta_k'*vec)
                      theta = 0.8*((delta_k'*vec)/((delta_k'*vec)- (gamma_kl'*delta_k)));
                  else
                      theta =1;
                  end
                  
              gamma_k = ((1-theta)*vec) + (theta*gamma_kl);
              k = k+1;
                  if k == 1
                      eta_1 = (norm(gamma_k)^2)/(gamma_k'*delta_k);
                      M = eta_1*eye(size(M));
                      
                      
                  else
                      M = M - ((vec*vec')/(delta_k'*vec)) + ((gamma_k*gamma_k')/(gamma_k'*delta_k));
                  end
                
                x = x_new;
                gf=gf_new;
                jc = jc_new;
                c=c_new;
                ci=ci_new;
                ai=ai_new;
                lme = lm_new;
                lmi=lmi_new;
                norm_gl = max(abs(gl_new));
                norm_ce = max(abs(c));
                norm_diff = max(abs(min(lmi_new,-ci)));
                condi = cond(M);
                deltaF = [deltaF max([norm_diff,norm_ce,norm_gl])];
                figure(2);
                clf;
                plot(1:k-1,log(deltaF(2:k)./deltaF(1:(k-1))));
                
              
              
              if options.verb == 1
                fprintf('%4d  %11.4e %11.4e %14d  %8.1e  %8.1e  %8.1e %8.1e %8.1e \n', ...
                k, max(abs(gl_new)), max(abs(c)),norm_diff, max(abs(lme)), max(abs(x)),alpha,theta,condi)
              end
            end
        else
           while k<= options.maxit

              if [norm_gl,norm_ce,norm_diff]<= options.tol
                 if options.verb == 1
                        fprintf('----------------------------------------------------------------------\n')
                 elseif options.verb == 2
                        fprintf('--------------------------------------------------------------\n')
                 end
                 break
              end
                        %sinon on résout (4.2)
              param = optimset('Display','off');
              [dk, fval, exitflag,output, lambda] = quadprog(M,gf,[],[],jc,-c,[],[],[],param);

              if (options.verb == 1) && (k == 0)
                    fprintf('\n')
                    fprintf('algorithme d optimisation quadratique successive sans contraintes d égalité: \n')
                    fprintf('----------------------------------------------\n')    
                    fprintf('%4s  %10s  %10s %10s %7s  %7s %7s \n',...
                            'iter', '   |gl|   ', '   |ce|   ','   ci<=0   ', '  |x|  ', ' |lme| ','   |lmi|   ')
              end
              %calcul des nouveaux x, lambda:
              x = x + dk;
              lm = lambda.eqlin;
              [f,c, ci, gf, jc,ai, hl, indic] = simul(6, x, lm,lmi);

              if strcmp(lmi,'no')==0
                  gl = gf + jc'*lm+ai'*lmi ;
                  [L, d] = cholmod(hl,small,big);
                  D = diag(d);
                  M = L*D*L';
                  size(M)
                  norm_gl = max(abs(gl));
                  norm_ce = max(abs(c));
                  norm_diff = max(abs(min(lmi,-ci)));
                  k = k+1;
              else
                  gl = gf + jc'*lm;
                  [L, d] = cholmod(hl,small,big);
                  D = diag(d);
                  M = L*D*L';
                  norm_gl = max(abs(gl));
                  norm_ce = max(abs(c));
                  norm_diff = 0;
                  k = k+1 ;
                  if options.verb == 1
                    fprintf('%4d  %10.4e %10.4e %10d %10.1e  %10.1e  %s\n', ...
                    k, max(abs(gl)), max(abs(c)), all(ci<=0), max(abs(x)), max(abs(lme)), 'NA')
                  end
              end
           end
        end        
    
    elseif options.deriv == 2 &&options.sqp == 1
            
        mi = length(lmi);
        small = 1.e-5;
        big = 1.e+5;
        
        % appel au simulateur
        [f,c,ci, gf, jc,ai, hl, indic] = simul(6, x, lme,lmi);
        
        if strcmp(lmi,'no')== 1
            gl = gf + jc'*lme 
        else 
            gl = gf + jc'*lme + ai'*lmi;
        end
        [L, d] = cholmod(hl,small,big);
        D = diag(d);
        M = L*D*L';
        tol = options.tol(1);
        norm_gl = max(abs(gl));
        norm_ce = max(abs(c));
        
        if strcmp(lmi,'no')==1
            norm_diff = 0;
        else    
            norm_diff = max(abs(min(lmi,-ci)));
        end
        if strcmp(lmi,'no')==0
            while k<= options.maxit

              if [norm_gl,norm_ce,norm_diff]<= options.tol
                  if options.verb == 1
                        fprintf('----------------------------------------------------------------------\n')
                  elseif options.verb == 2
                        fprintf('--------------------------------------------------------------\n')
                  end
                  break
              end

              %sinon on résout (4.2)
              param = optimset('Display', 'off');
              [dk, fval, exitflag,output, lambda] = quadprog(M,gf,ai,-ci,jc,-c,[],[],[], param);
              
              if (options.verb == 1) && (k == 0)
                fprintf('\n')
                fprintf('algorithme d optimisation quadratique successive: \n')
                fprintf('----------------------------------------------\n')    
                fprintf('%4s  %10s  %10s %11s %7s  %7s %7s \n',...
                        'iter', '   |gl|   ', '   |ce|   ','   ci<=0   ', '  |x|  ', ' |lme| ','   |lmi|   ')
              end
              %calcul des nouveaux x, lambda:
              x = x + dk;
              lm = lambda.eqlin;
              lmi = lambda.ineqlin;

              [f,c, ci, gf, jc,ai, hl, indic] = simul(6, x, lm,lmi);
              gl = gf + jc'*lm + ai'*lmi;
              [L, d] = cholmod(hl,small,big);
              D = diag(d);
              M = L*D*L';
              norm_gl = max(abs(gl));
              norm_ce = max(abs(c));
              norm_diff = max(abs(min(lmi,-ci)));
              k = k+1;
              if options.verb == 1
                fprintf('%4d  %10.4e %10.4e %10d  %7.1e  %7.1e  %7.1e\n', ...
                k, max(abs(gl)), max(abs(c)),all(ci<=0), max(abs(x)), max(abs(lme)),max(abs(lmi)))
              end
            end
        else
           while k<= options.maxit

              if [norm_gl,norm_ce,norm_diff]<= options.tol
                 if options.verb == 1
                        fprintf('----------------------------------------------------------------------\n')
                 elseif options.verb == 2
                        fprintf('--------------------------------------------------------------\n')
                 end
                 break
              end
                        %sinon on résout (4.2)
              param = optimset('Display','off');
              [dk, fval, exitflag,output, lambda] = quadprog(M,gf,[],[],jc,-c,[],[],[],param);

              if (options.verb == 1) && (k == 0)
                    fprintf('\n')
                    fprintf('algorithme d optimisation quadratique successive sans contraintes d égalité: \n')
                    fprintf('----------------------------------------------\n')    
                    fprintf('%4s  %10s  %10s %11s %7s  %7s %7s \n',...
                            'iter', '   |gl|   ', '   |ce|   ','   ci<=0   ', '  |x|  ', ' |lme| ','   |lmi|   ')
              end
              %calcul des nouveaux x, lambda:
              x = x + dk;
              lm = lambda.eqlin;
              [f,c, ci, gf, jc,ai, hl, indic] = simul(6, x, lm,lmi);

              if strcmp(lmi,'no')==0
                  gl = gf + jc'*lm+ai'*lmi ;
                  [L, d] = cholmod(hl,small,big);
                  D = diag(d);
                  M = L*D*L';
                  size(M)
                  norm_gl = max(abs(gl));
                  norm_ce = max(abs(c));
                  norm_diff = max(abs(min(lmi,-ci)));
                  k = k+1
              else
                  gl = gf + jc'*lm;
                  [L, d] = cholmod(hl,small,big);
                  D = diag(d);
                  M = L*D*L';
                  norm_gl = max(abs(gl));
                  norm_ce = max(abs(c));
                  norm_diff = 0;
                  k = k+1 ;
                  if options.verb == 1
                    fprintf('%4d  %10.4e %10.4e %10d %7.1e  %7.1e  %s\n', ...
                    k, max(abs(gl)), max(abs(c)), all(ci<=0), max(abs(x)), max(abs(lme)), 'NA')
                  end
              end
           end
        end
       
        
        
    elseif options.deriv == 2 && options.sqp == 2
    % boucle de l'algo de Newton avec recherche lineaire
    while ((k < options.maxit) && ...
           (indic == 0)        && ...
           (options.rl == 0)   )
        
        % appel au simulateur :
        %   - f   fonction a minimiser
        %   - ce  contraintes d'egalite
        %   - gf  gradient de f
        %   - jc  jacobienne des contraintes
        %   - hl  hessienne du Lagrangien
        [f, ce,ci, gf, jc,ai, hl, indic] = simul(6, x, lme,lmi);
        nsim = nsim + 1;
        
        gl = gf + jc'*lme;    % gradient du Lagrangien
        
        % test de la condition d'optimalite
        if (max(abs(gl)) <= options.tol(1)) || (max(abs(ce)) <= options.tol(2))
          
            flag = 1;
            
            % affichage
            if options.verb == 1
                fprintf('----------------------------------------------------------------------\n')
            elseif options.verb == 2
                fprintf('--------------------------------------------------------------\n')
            end
            
            break
        end
        
        % resolution lineaire de l'iteration de Newton
        F = [gl ; ce];
        jF = [[hl jc'];[jc zeros(m)]];
        p = - full(jF)\F;
        
        %parametres pour la recherche lineaire
        phi = (1./2)*F'*F;
        pente = -F'*F;
        alpha = 1.;     % pas d'Armijo
        i = 0;          % compteur d'iterations d'Armijo
        
        % affichage
        if options.verb == 1 && k == 0
            fprintf('\n')
            fprintf('algorithme de Newton avec recherche lineaire: \n')
            fprintf('----------------------------------------------------------------------\n')    
            fprintf('%4s  %10s  %10s  %7s  %7s  %9s  %11s \n',...
                    'iter', '   |gl|   ', '   |ce|   ', '  |x|  ', ' |lme| ', '  alpha  ', '    phi    ')
        end
        
        % affichage
        if options.verb == 2
            if k == 0
                fprintf('\n')
                fprintf('algorithme de Newton avec recherche lineaire: \n')
            end
            fprintf('--------------------------------------------------------------\n')
            fprintf('iter %3d,  simul %3d,  phi %12.5e,  pente %12.5e \n', ...
                    k+1, nsim, phi, pente)
            fprintf('\n')
            fprintf('  recherche lineaire d Armijo'': |d| = %9.2e \n', max(abs(p(n+1:n+m))) )
            fprintf('  %10s  %12s  %12s \n',...
                    '  alpha   ', '  phip-phi  ', '   DF(phi)  ')
        end
        
        % boucle de la recherche lineaire avec regle d'Armijo
        while 1
            
            % test de depassement du nombre max d'iterations
            if i > options.maxit
                fprintf('\n')
                fprintf('error Armijo : no valid alpha found \n')
                indic = 1;
                break
            end
            
            % nouvelles variables a tester
            new_x = x + alpha * p(1:n);
            new_lme = lme + alpha * p(n+1:n+m);
            
            [new_f, new_ce,new_ci, new_gf, new_jc,new_ai, new_hl, indic] = simul(6, new_x, new_lme,'no');
            nsim = nsim + 1;
            
            % test d'erreur du simulateur
            if indic == 1
                break
            end
            
            new_gl = new_gf + new_jc'*new_lme;
            new_F = [new_gl ; new_ce];
            
            phip = (1./2)*new_F'*new_F;
            dphi = phip - phi;
            Dphi = dphi / alpha;
            
            % affichage
            if options.verb == 2
                fprintf('  %10.4e  %12.5e  %12.5e \n',...
                        alpha, dphi, Dphi)
            end
            
            % test de la regle d'Armijo
            if Dphi <= omega * pente
                
                % l'iteration est validee, on peut incrementer le compteur
                k = k+1;
                
                % affichage
                if options.verb == 1
                    fprintf('%4d  %10.4e  %10.4e  %7.1e  %7.1e  %9.3e  %11.5e \n', ...
                            k, max(abs(gl)), max(abs(ce)), max(abs(x)), max(abs(lme)), alpha, phi)
                elseif options.verb == 2
                    fprintf('\n')
                    fprintf('  |gl| = %9.3e,  |ce| = %9.3e \n', max(abs(gl)), max(abs(ce)))
                    fprintf('\n')
                end
                
                % mise a jour des variables pour la prochaine etape de boucle
                x = new_x;
                lme = new_lme;
                deltaF = [ deltaF [ max(abs(new_F - F)) ] ];
                
                break
            end
            
            % mise a jour des parametres d'Armijo
            alpha = alpha / 2.;
            i = i+1;
            
        end
        
    end
    
    % boucle de l'algo de Newton sans recherche lineaire
    while ((k < options.maxit) && ...
           (indic == 0)        && ...
           (options.rl == 1)   )
        
        % appel au simulateur
        [f, ce,ci, gf, jc,ai, hl, indic] = simul(6, x, lme,lmi);
        nsim = nsim + 1;
        
        % gradient du Lagrangien
        gl = gf + jc'*lme;
        
        % test de la condition d'optimalite
        if (max(abs(gl)) <= options.tol(1)) || (max(abs(ce)) <= options.tol(2))
          
            flag = 1;
            
            % affichage
            if options.verb == 1
                fprintf('----------------------------------------------\n')
            end
            
            break
        end
        
        % resolution lineaire de l'iteration de Newton
        F = [gl ; ce];
        jF = [[hl jc'];[jc zeros(m)]];
        p = -jF\F;
        
        % affichage
        if options.verb == 2
            fprintf('\n')
            fprintf('sqp warning : pas de recherche lineaire ici \n')
            options.verb = 1;
        end
        
        % affichage
        if options.verb == 1 && k == 0
            fprintf('\n')
            fprintf('algorithme de Newton sans recherche lineaire: \n')
            fprintf('----------------------------------------------\n')    
            fprintf('%4s  %10s  %10s  %7s  %7s \n',...
                    'iter', '   |gl|   ', '   |ce|   ', '  |x|  ', ' |lme| ')
        end
        
        % variables pour la prochaine etape
        new_x = x + p(1:n);
        new_lme = lme + p(n+1:n+m);
            
        [new_f, new_ce,new_ci, new_gf, new_jc,new_ai, new_hl, indic] = simul(6, new_x, new_lme,'no');
        nsim = nsim + 1;
            
        new_gl = new_gf + new_jc'*new_lme;
        new_F = [new_gl ; new_ce];
        
        %mise a jour du compteur        
        k = k+1;
        
        % affichage
        if options.verb == 1
            fprintf('%4d  %10.4e  %10.4e  %7.1e  %7.1e  \n', ...
                    k, max(abs(gl)), max(abs(ce)), max(abs(x)), max(abs(lme)))
        end
        
        % mise a jour des variables       
        x = new_x;
        lme = new_lme;
        deltaF = [ deltaF [ max(abs(new_F - F)) ] ];
        
    end
    
    
    % affichage graphique de l'etat final de la recherche
    [f, ce,ci, gf, jc,ai, hl, indic] = simul(2, x, lme,lmi);
      
    figure(1);
    grid on;
    shg;
    figure(2);
    clf;
    plot(1:k,deltaF(1:k));
    figure(3);
    clf;
    plot(floor(k/2):k,deltaF(floor(k/2):k));
    end
    % sauvegarde de l'etat final de l'algorithme et affichage des erreurs
    if flag
        info.status = 0;
      
    elseif indic == 1
        fprintf('\n')
        fprintf('error chs : inconsistant input \n')
        info.status = 1;
      
    elseif k >= options.maxit
        fprintf('\n')
        fprintf('warning spq : depassement options.maxit \n')
        info.status = 2;
    end
    info.niter = k;
    
end