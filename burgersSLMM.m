%function burgersSLMM()
% Solving Burgers with a SL method with moving meshes.
%
% This uses the following external files.
% ./mmpde5.m
% ./fwd_euler.m
% ~/dos/MATLAB/Eqd1dExact.m
% ~/dos/MATLAB/interp1cubicL.m
% ~/dos/MATLAB/ppval_lim.m

clear all
clf

% Parameters
% General Parameters
N=21;       % Level of spacial discretisation
K=40;         % Departure point iterations
%Dt= 0.012;   % Timestep
tN = 150;
tmax = 1.5;   % Final time
theta_t=1/2;  % Theta for the Theta method in time.
theta_x=1/2;  % Theta for the Theta-method for departure points.
epsilon=0.0001 ; % Epsilon in the PDE

Dt = tmax/tN;

% Mesh Parameters
%mesh = 'static';
%mesh = 'prescribed';
mesh = 'moving-exact';  % Mesh movement type
%mesh = 'moving-relax';  % Mesh movement type
limiter = 1;        % Flux limiter for interpolation
interpolation = 'linear';
%interpolation = 'CSpline';
%interpolation = 'CLagrange';
%interpolation = 'pchip';

b = 1;                 % Parameter for the monitor function.
m=@(x,u,uprime) sqrt(b + uprime.^2);
%m=@(x,u,uprime)ones(size(u))
p_smooth = 5;
tau =1;
with_euler = 1;



% Initial and boundary conditions
u0 = @(x) (sin(2*pi*x) + 1/2*sin(pi*x));
%u0 = @(x) 1 - tanh(x/(2*epsilon));
%u_l = 1;
u_l = 0;
u_r = 0;
Dx= 1./(N+1); % Space step

% Matrices to be used. Tridiagonal, so could be done more efficiently

% del2 - \delta_x^2, a tridiagonal matrix, finite difference second 
% derivative operator.
del2 = 2*eye(N) - diag(ones(N-1,1),1) - diag(ones(N-1,1),-1);
% LHS and RHS of SISL formulation
M_RHS = eye(N) - Dt/(Dx^2) * (1 - theta_t) * epsilon * del2;
M_LHS = eye(N) + Dt/(Dx^2) * theta_t * epsilon * del2;

% Initialisation
X0 = 0;
X = (1:N)'./(N+1);
XN = 1;
Dxi = 1/(N+1);
%X = X.^2;
Un = u0(X);
X_An = X;
TT = (0:Dt:tmax)';
XX = zeros(length(TT),length(X));
XX(1,:) = X;
dept_convergence = zeros(K,length(TT));
uout= XX;
uout(1,:) = Un;
%jj=1; % Counting variable for saving entries to uout.

del2n1 = eye(N);
del2n = eye(N);

% Outer loop. Timestep
for tt = 1:length(TT)
    t = TT(tt);
% Timestep initialisation
% Cubic spline of Un (including the constant end points)  

% Mesh movement
switch mesh
    case 'static'
        X_An1 = X;
	vtitle = ['Semi-Lagrangian Burgers, static mesh, N = '...
	,num2str(N)];
    case 'prescribed'
        % prescribed mesh, loaded from file, innit.
	% Not yet done, static for now.
        X_An1 = X;
    case 'moving-exact'
        %XN = [0;X_An;1];
        U = [u_l;Un;u_r];
        DUDX = [u_l;diff(U)]./[1;diff([0;X_An;1])];
        if with_euler
            Uhat = U + Dt*fwd_euler(U,[0;X_An;1],epsilon);
            DUhatDX = [u_l;diff(Uhat)]./[1;diff([0;X_An;1])];
            M = m([0;X_An;1],Uhat,DUhatDX);
        else
            M = m([0;X_An;1],U,DUDX);
        end
        % Smooth
        for iii = 1:p_smooth
            M = [2/3*M(2)+1/3*M(1); M; 2/3*M(end-1) + 1/3*M(end)];
            M = conv(M,[1/4,1/2,1/4],'same');  % Smoothing
            M = M(2:end-1);
        end % iii
        M=M./(trapz([0;X_An;1],M)); % normalise
        M = 1/2 + M./2;             % Average
        %
        monitor.type = 'points';
        monitor.x = [0;X_An;1];
        monitor.M = M;
        X_An1 = Eqd1dExact([0;X_An;1], monitor);
        X_An1 = X_An1(2:end-1);
 %       X_An1 = (X_An1 + X_An)./2;
	vtitle = ['Semi-Lagrangian Burgers, moving mesh (exact), N = '...
	,num2str(N)];
    case 'moving-relax'
        % This is different from the above method, in that 
        % we solve a moving mesh PDE to get the movement, 
        % and to do the we take M to sit at the midpoints.
        
        %XN = [0;X_An;1];
        U = [u_l;Un;u_r];
        DUDX = diff(U)./diff([0;X_An;1]);    % ! Sits at the midpoints
        M = m([0;X_An;1],U,DUDX); % This might require changing for
                                  % different monitor functions.
        % Smooth
        for iii = 1:p_smooth
            M = [2/3*M(2)+1/3*M(1); M; 2/3*M(end-1) + 1/3*M(end)];
            % This line is needed to correct the smoothing stencil
            % at the endpoints ( change from [1/4 1/2 1/4] which will
            % use the zero, to [1/3 2/3 0] ).
            M = conv(M,[1/4,1/2,1/4],'same');  % Using convolution
            M = M(2:end-1); % Undo the first line of this loop.
        end % iii
        M=M./(trapz(1/2*([0;X_An]+[X_An;1]),M)); % normalise
        M = 1/2 + M./2;             % Average
        %
        x = [0;X_An;1];
        %monitor.M = m([0;X_An;1],U,DUDX);
        % Spline it!
        M_pp = spline(x,mmpde5(M,x,Dxi,tau));
        [~,ode_out] = ode15s(@(t,x_)ppval(M_pp,x_),[0 Dt],x);
        X_An1 = ode_out(end,:)';
        %X_An1 = [0;X_An;1] + Dt*mmpde5(M,x,Dxi,tau);
        X_An1 = X_An1(2:end-1);
 %       X_An1 = (X_An1 + X_An)./2;
	vtitle = ['Semi-Lagrangian Burgers, moving mesh (relax), N = '...
	,num2str(N)];
    otherwise
        disp('Unknown movement type, using static')
        X_An1 = X;
end % switch

DX_An = diff([0;X_An;1]);
DX_An1 = diff([0;X_An1;1]);

del2n = del2n1;

% Build up the LHS and RHS tridiagonal matrix. Can recycle these 
% and get 2 uses out of each del2 matrix (only redo one)
del2n1(1,1) = -2/(DX_An1(2)*DX_An1(1));
del2n1(1,2) = 2/(DX_An1(2)*(DX_An1(2)+DX_An1(1)));
for ii=2:(N-1)
    del2n1(ii,ii-1) = 2/(DX_An1(ii)*(DX_An1(ii+1)+DX_An1(ii)));
    del2n1(ii,ii) = -2/(DX_An1(ii)*DX_An1(ii+1));
    del2n1(ii,ii+1) = 2/(DX_An1(ii+1)*(DX_An1(ii+1)+DX_An1(ii)));
end
del2n1(N,N-1) = 2/(DX_An1(N-1)*(DX_An1(N)+DX_An1(N-1)));
del2n1(N,N) = -2/(DX_An1(N)*DX_An1(N-1));
M_LHS = eye(N) - Dt * theta_t * epsilon * del2n1;

%del2n(1,1) = -2/(DX_An(2)*DX_An(1));
%del2n(1,2) = 2/(DX_An(2)*(DX_An(2)+DX_An(1)));
%for ii=2:(N-1)
%    del2n(ii,ii-1) = 2/(DX_An(ii)*(DX_An(ii+1)+DX_An(ii)));
%    del2n(ii,ii) = -2/(DX_An(ii)*DX_An(ii+1));
%    del2n(ii,ii+1) = 2/(DX_An(ii+1)*(DX_An(ii+1)+DX_An(ii)));
%end
%del2n(N,N-1) = 2/(DX_An(N-1)*(DX_An(N)+DX_An(N-1)));
%del2n(N,N) = -2/(DX_An(N)*DX_An(N-1));
M_RHS = eye(N) + Dt * (1 - theta_t) * epsilon * del2n;

% Make the spline for interpolation
rhs0 = -Dt/(Dx^2) * (1-theta_t) * epsilon *...
    (5*Un(1) + 4*Un(2) - Un(3)); % This is an order h^2 approx'n to u''
rhsN = -Dt/(Dx^2) * (1-theta_t) * epsilon *...
    (5*Un(N) + 4*Un(N-1) - Un(N-2));
switch interpolation
    case 'linear'
        pp_rhs = interp1([0;X_An;1],[rhs0; (M_RHS * Un) ; rhsN],...
            'linear','pp');
    case 'CSpline'
        pp_rhs = spline([0;X_An;1],[rhs0; (M_RHS * Un) ; rhsN]);
    case 'CLagrange'
        pp_rhs = interp1cubicL([0;X_An;1],[rhs0; (M_RHS * Un) ; rhsN]);
    case 'pchip'
        pp_rhs = pchip([0;X_An;1],[rhs0; (M_RHS * Un) ; rhsN]);
end

% Initial guess of departure points.
X_D = X_An1 - Dt*Un;
 % Inner loop. 
 for k = 1:K
 % Evaluate RHS at the departure points.
 if limiter
     rhs_D = ppval_lim(pp_rhs,X_D);
 else
     rhs_D = ppval(pp_rhs, X_D);
 end
 
 % Solve the implicit equation (Thomas algorithm implemented later)
 U_A = M_LHS\rhs_D;
 % New guess at the departure points (theta-method) and recalculate 
 % RHS_D and U_A. 
 X_D_old = X_D;
 X_D = X_An1 - Dt*(theta_x*U_A + (1-theta_x)*Un);
 dept_convergence(k,tt) = norm(X_D - X_D_old);
 %rhs_D = ppval(pp_rhs, X_D);
 if limiter
     rhs_D = ppval_lim(pp_rhs,X_D);
 else
     rhs_D = ppval(pp_rhs, X_D);
 end % if limiter
 end % for k
 U_A = M_LHS\rhs_D;
 % End of Inner Loop


%%% Prep for next timestep. %%%
X_An = X_An1;

Un = U_A;
%jj=jj+1;
uout(tt,:) = Un;
%%% And plot %%% 
%plot([0;X_An1;1],[u_l;U_A;u_r]), title(['t = ',num2str(t)]), ylim([-1 1.5])
%if isequal(mesh,'moving')
%    hold on
%    plot(monitor.x,monitor.M,'g-')
%    hold off
%end
%drawnow()
XX(tt,:) = X_An;
bigX_D(tt,:) = X_D;
end % for t
%pause
for i = 1:N
    plot(XX(:,i),TT)
    hold on
    %plot(bigX_D(:,i),TT,'g-')
end

% Okay, only want ~200 frames in the movie.
t_skip = ceil(length(TT)/200);
mk_video('test.avi',vtitle,TT(1:t_skip:end),uout(1:t_skip:end,:),...
    XX(1:t_skip:end,:))

%end % function main

