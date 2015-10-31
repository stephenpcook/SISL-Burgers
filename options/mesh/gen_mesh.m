clear all

%% moving_exact_arc
mesh_movement = 'moving-exact';  % Mesh movement type
b = 0.1;
m=@(x,u,u_x,u_xx) sqrt(b + u_x.^2);
p_smooth = 5;
with_euler = 0;

save('moving_exact_arc'...
  ,'mesh_movement'...
  ,'m', 'p_smooth', 'with_euler'...
  );

%% moving_exact_curv
mesh_movement = 'moving-exact';  % Mesh movement type
b = 0.1;
m=@(x,u,u_x,u_xx) sqrt(b + u_xx.^2);
p_smooth = 5;
with_euler = 0;

save('moving_exact_curv'...
  ,'mesh_movement'...
  ,'m', 'p_smooth', 'with_euler'...
  );
