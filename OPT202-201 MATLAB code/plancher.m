% Dessin du plancher

  global R S    % les p fonctions affines décrivant le plancher "max(R+S*x)"
  global FIG_XMIN FIG_XMAX FIG_YMIN FIG_YMAX	% limites de la figure

  floor_color = 0.81*[1 1 1];	% couleur du plancher

  % dessin du plancher

  xf = [FIG_XMIN;FIG_XMIN;FIG_XMAX;FIG_XMAX];
  for i=1:length(R)
    yf = [R(i)+S(i)*FIG_XMIN;FIG_YMIN;FIG_YMIN;R(i)+S(i)*FIG_XMAX];
    fill(xf,yf,floor_color,'EdgeColor',floor_color);
  end
