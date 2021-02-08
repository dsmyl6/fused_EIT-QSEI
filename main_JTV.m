%%%%% Sample 2D JTV code (ERT + Elasticity Imaging) %%%%%
%%% Coded by: Danny Smyl 
%%% Note: hyperparameters need adjustment depending on discretization 
%%% Disclaimer 1: Code is not optimized and superfluous variables may exist
%%% Disclaimer 2: EIDORS is needed to execute ERT functions 
%%% Disclaimer 3: Inversion mesh is very coarse, finer mesh = less model error
clear all
addpath (genpath ('C:\Users\danny\Desktop\JTV')) %include filepath
run C:\Users\danny\Dropbox\Danny\AALTO\ERT\ERT_Base\eidors-v3.9.1-ng\eidors\startup %intialize EIDORS
% Add EIDORS filepath if needed

%%%%%%%%%%%%%%%%%%%%%%%
%%% Mesh parameters %%%
%%%%%%%%%%%%%%%%%%%%%%%
squarewidth  = 10;
squareheight = 10;
L_x          = 4;% #Electrodes horizontally
L_y          = 4;% #Electrodes vertically
el_width     = 1;% Electrode width
el_depth     = 1;% Electrode depth
L = 2*(L_x + L_y);% Number of electrodes

%%% Mesh Params %%%
xy = [squarewidth squareheight; 0 squareheight; 0 0; squarewidth 0]; %outer boundary of mesh
es = [el_width el_depth];

el_pos1= [[squarewidth/5 2*squarewidth/5 3*squarewidth/5 4*squarewidth/5]'   zeros(L_y,1)];
el_pos2= [squarewidth*ones(L_x,1) [squarewidth/5 2*squarewidth/5 3*squarewidth/5 4*squarewidth/5]'];
el_pos3= [[4*squarewidth/5 3*squarewidth/5 2*squarewidth/5 squarewidth/5]' , squareheight*ones(L_y,1)];
el_pos4= [zeros(L_x,1), [4*squarewidth/5 3*squarewidth/5 2*squarewidth/5 squarewidth/5]'];
el_pos = [el_pos1;el_pos2;el_pos3;el_pos4];

fmdlF = ng_mk_2d_model({xy, 0.2*el_width}, el_pos, es);
show_fem(fmdlF,[0 1]),axis off;

%%% Fine (simulation data) mesh
fmdlF.solve =  @fwd_solve_higher_order;
fmdlF.system_mat = @system_mat_higher_order;
fmdlF.jacobian = @jacobian_adjoint_higher_order;
fmdlF.approx_type    = 'tri6';

gF = fmdlF.nodes; %nodal coords
HF = fmdlF.elems;
sNF = length(HF(:,1));%size of param space

point = [];
for q = 1:sNF
    gx(q,:) = gF(HF(q,:),1);
    gy(q,:) = gF(HF(q,:),2);
    
    gcx(q) = mean(gx(q,:));
    gcy(q) = mean(gy(q,:));
end
gF = [gcx; gcy]'; %element center coords

%%% plot ERT simulation mesh%%%
figure(1)
show_fem(fmdlF); axis off
drawnow

%%% ERT Stimulation pattern %%%
I = toeplitz([1;-1;zeros(L-2,1)],...
    [1,zeros(1,L-2),-1]);
I = 2*I/max(max(I)); %change to 2mV
MeasPatt = toeplitz([1;-1;zeros(L-2,1)],[1,zeros(1,L-2),-1]);

for jj = 1:L
    stim(jj).stim_pattern = I(:,jj);
    stim(jj).meas_pattern = MeasPatt';
end

fmdlF.stimulation = stim;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Generate FEM mesh for QSEI %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Lx = squarewidth;
Ly = squareheight;
elwidth = Ly/40;

%%% Elasticity simulation mesh %%%
[gQSEI,xQSEI,yQSEI,ExQSEI,EyQSEI,TriQSEI,constraint,nelQSEI,nelxQSEI,nelyQSEI] = meshgen(Lx,Ly,elwidth);

%%% ERT data generation mesh%%%
g0 = gF; % nodes
H0 = HF; % connectivity

%%% Elastic params (if random crack is not generated %%%
Emax = 200;
Emin = 10;
nu = 0.3; %poisson
th = 0.02; %thickness

%%% conductivity params %%%
sigmax = 20;

%%% conductivity defect generation %%%
select_fcn = inline('(x-2.5).^2+(y-5).^2<0.6^2','x','y','z');
im_target = mk_image(fmdlF,sigmax);
im_target.elem_data = sigmax - 0.99999*sigmax*elem_select(im_target.fwd_model, select_fcn);

%%% elastic defect
cE = knnsearch([ExQSEI,EyQSEI],[2.5,5],'K',60); %elastic defect

Ehomo = 200; %homgeneous elastic background from prior material info
Ehomogeneous = Ehomo*ones(length(ExQSEI),1);
Etrue = 200*ones(length(ExQSEI),1);
Etrue(cE) = Emin; % set crack near 0

%%% Plotting target conductivity field %%%
IO = isInterior(TriQSEI);

figure(2),
show_fem(im_target),title('$\sigma_{true}$','FontSize',14); %target damage image

%%%%%%%%%% Plotting target elastic field%%%%%%%%%%
figure(3), 
F = scatteredInterpolant(ExQSEI,EyQSEI,Etrue);
int = F(gQSEI(:,1),gQSEI(:,2));
trisurf(TriQSEI(IO,:),TriQSEI.Points(:,1), TriQSEI.Points(:,2),int),view(2),colormap('hot'),set(0,'defaulttextInterpreter','latex'),daspect([1 1 1]), colorbar('eastoutside');
axis tight,title('$E_{true}$','FontSize',14), axis off, box on;
drawnow

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% SIMULATE DATA %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
UF = fwd_solve(im_target);
Uel2_nonoise = UF.meas; %noiseless ERT data

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Simulate QSEI measurements %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  boundary conditions for Elastic inv model %%%
%%% Enternal forces (end)%%%
Fend0 = 1e-1; %%% Total (non-distributed) end load
force1 = find(gQSEI(:,1) > L_x - 10e-6);
Fend = Fend0/numel(force1); %distributed load on right side
force = [force1 ones(length(force1),1) Fend*ones(length(force1),1)];

%%% support conditions (fixed left side, restrained right side (y-dir)
cond1 = [find(gQSEI(:,1) < 10e-6),find(gQSEI(:,1) < 10e-6)];
cond2 = [ones(1,length(cond1)); 2*ones(1,length(cond1))];
cond3 = [force1 2*ones(length(force1),1)];
cond = [cond1(:) cond2(:)]; cond = [cond; cond3]; condF = cond;

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Sim data (Elastic) %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
[u,ux,uy,K,Ael,~] = SimulateData(Etrue,nu,th,constraint,cond,force,gQSEI(:,1),gQSEI(:,2),TriQSEI);
um = u; %noiseless QSEI data

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% coarse mesh load for inverse problem %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% ERT Coarse  Inversion mesh%%%%%%%
fmdlC = ng_mk_2d_model({xy, 0.5*el_width}, el_pos, es);

figure(4)
show_fem(fmdlC); axis off
drawnow

fmdlC.solve =  @fwd_solve_higher_order;
fmdlC.system_mat = @system_mat_higher_order;
fmdlC.jacobian = @jacobian_adjoint_higher_order;
fmdlC.approx_type    = 'tri6';

gC = fmdlC.nodes;
HC = fmdlC.elems;
sN = length(HC(:,1));

point = [];
for q = 1:sN
    gx(q,:) = gC(HC(q,:),1);
    gy(q,:) = gC(HC(q,:),2);
    
    gcx(q) = mean(gx(q,:));
    gcy(q) = mean(gy(q,:));
end
g = [gcx; gcy]';

%%% set stim pattern %%%
fmdlC.stimulation = stim;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Measurement noise ERT %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%
meas_noise_coef = 1e-3;  %Magnitude of measurement noise

%%% add noise to EIT measurements %%%
Uel = Uel2_nonoise + meas_noise_coef*abs(Uel2_nonoise).*randn(size(Uel2_nonoise));
Uel = Uel + (meas_noise_coef*(max(Uel)-min(Uel)))*randn(size(Uel));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Remesh QSEI (elastic inverse mesh %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
TriQSEIinv = delaunayTriangulation(gC(:,1),gC(:,2)); %Tri is H
Eloc = incenter(TriQSEIinv);
ExQSEIinv= Eloc(:,1);
EyQSEIinv= Eloc(:,2);
nelQSEIinv=length(EyQSEIinv);
gQSEIinv = gC;
IO = isInterior(TriQSEIinv); %for plotting

%%%  boundary conditions and forces %%%
force1 = find(gQSEIinv(:,1) > L_x - 10e-6);
Fend = Fend0/numel(force1); %distributed load on right side
force = [force1 ones(length(force1),1) Fend*ones(length(force1),1)];

%%% support conditions (fixed left side, restrained right side (y-dir)
cond1 = [find(gQSEIinv(:,1) < 10e-6),find(gQSEIinv(:,1) < 10e-6)];
cond2 = [ones(1,length(cond1)); 2*ones(1,length(cond1))];
cond3 = [force1 2*ones(length(force1),1)];
cond = [cond1(:) cond2(:)]; cond = [cond; cond3]; condF = cond;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Displacement Mapping (interp) from Fine to Coarse grid  %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% DAMAGED DATA interpolation onto coarse inversion mesh %%%
F=scatteredInterpolant(gQSEI(:,1),gQSEI(:,2),ux','natural');
uxinv = F(gQSEIinv(:,1),gQSEIinv(:,2));

F=scatteredInterpolant(gQSEI(:,1),gQSEI(:,2),uy','natural');
uyinv = F(gQSEIinv(:,1),gQSEIinv(:,2));

%%% odd are x displacements
um = [];

for gg = 1:(length(uxinv))
    if gg == 1
        um = [uxinv(gg); uyinv(gg)];
    else
        umd = [uxinv(gg); uyinv(gg)];
        um = [um; umd];
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Add Measurement noise QSEI %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
um = um + meas_noise_coef*abs(um).*randn(size(um));
um = um + (meas_noise_coef*(max(um)-min(um)))*randn(size(um));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% homogenous estimate for the E %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Eexp = Emax*ones(numel(ExQSEIinv),1); % Expected value for E
sigmaexp = sigmax*ones(sN,1); % Expected value for sigma
sig = sigmaexp; %starting point
E_GN = Eexp; %starting point

%%%%%%%%%%%%%%%%%%%%
%%%% JTV Prior  %%%% 
%%%%%%%%%%%%%%%%%%%%
alpha = 1e-9; %Needs adjusting depending on discretization
beta = 0.0001; %Needs adjusting depending on discretization
alpha_E = sigmaexp(1)/Eexp(1);

%%% Tik Regularization matrix %%%
lambda1 = 1e-10;
lambda2 = lambda1; %if variable L_1 regularization is desired
RegLambda = lambda1*eye(sN);
RegLambda1 = sparse(RegLambda'*RegLambda);
RegLambda2 = sparse(RegLambda'*RegLambda);

%%% TV Regularization matrix %%%
imgGN = mk_image(fmdlC, sig);
L_TV= prior_TV(imgGN);

%%%%%%%%%%%%%% opt params %%%%%%%%%%%%%%%%%%
theta = [sig;E_GN];%stacked opt parameter
minval=10^-8; %proj val
ii = 1;

%%% initialize ERT FMDL %%%
imgGN = mk_image(fmdlC, sig(:,ii));
Uref = fwd_solve(imgGN);
Urefel = Uref.meas;

%%% initialize ERT FMDL %%%
[usim,~] = FMDL(E_GN(:,ii),nu,th,constraint,TriQSEIinv,gQSEIinv,cond,force);

%%% Gauss-Newton iterations %%%
Niter = 10;

for ii = 2:Niter
    
    %%%% ERT Jacobian (Adjoint) %%%%
    J_ERT = calc_jacobian(fmdlC,imgGN);
    
    %%%% QSEI Jacobian (complex differentiation) %%%%
    J_QSEI = PertubedJ(sN,E_GN(:,ii-1),nu,th,constraint,TriQSEIinv,gQSEIinv,cond,force);
       
    %%% JTV Prior %%%
    phi = sqrt( (L_TV*sig(:,ii-1)).^2 + (L_TV*(alpha_E*E_GN(:,ii-1)) ).^2 + beta);
    phi_inv = spdiags( 1./phi, 0, length(phi), length(phi));
    
    grad_TV_ERT = alpha*L_TV'*phi_inv*L_TV*(sig(:,ii-1));
    Hess_TV_ERT = alpha*L_TV'*phi_inv*L_TV;
    grad_TV_QSEI= alpha*L_TV'*phi_inv*L_TV*(alpha_E*E_GN(:,ii-1));
    Hess_TV_QSEI= Hess_TV_ERT;
    
    %%% ERT minimizer %%%
    HH1 = J_ERT'*J_ERT + Hess_TV_ERT + RegLambda1;
    zz1 =  J_ERT'*(Uel-Urefel) - grad_TV_ERT;
    dtheta_ERT = HH1\zz1;

    %%% QSEI minimizer %%%
    HH2 = J_QSEI'*J_QSEI + Hess_TV_QSEI + RegLambda2;
    zz2 =  J_QSEI'*(um - usim) - grad_TV_QSEI;
    dtheta_QSEI = HH2\zz2;
    
    %%% GN search direction [stacked] %%%
    dtheta = [dtheta_ERT; dtheta_QSEI];
    
    %%% basic linesearch %%%
    clear sk costfunct;
    for q = 1:41 % k \in (0,2]
        if q==1
            sk = 10^-7;
        else
            sk(q) = 0.05*(q-1);%0.05 is the linesearch step size
        end
        %%% ERT sim %%%
        sigsk = sig(:,ii-1) + sk(q)*dtheta_ERT;
        imgGNsk = mk_image(fmdlC, sigsk);
        Urefsk = fwd_solve(imgGNsk);
        Urefelsk = Urefsk.meas;
        %%% QSEI sim %%%
        Esk = E_GN(:,ii-1)+sk(q)*dtheta_QSEI;
        [usk,~] = FMDL(Esk,nu,th,constraint,TriQSEIinv,gQSEIinv,cond,force);
        %%% linesearch cost function %%%
        costfunct(q) = norm(um - usk)^2 + norm(Uel - Urefelsk)^2 + alpha*sum(sqrt((L_TV*(alpha_E*Esk)).^2 + (L_TV*sigsk).^2)) + norm(RegLambda2*(Esk),1) + norm(RegLambda1*(sigsk),1);
        
        figure(8);
        plot(costfunct);title('linesearch');
        drawnow
        
        if q > 1 && costfunct(q) > costfunct(q-1)
            break;
        end
    end
    
    [~,idxsk] = min(costfunct); % find optimal [minimal] step size
    k = sk(idxsk); % linesearch stepsize
    
    %%%%%%%%%%%%%%%%%% GN Update %%%%%%%%%%%%%%%%%%%%
    theta(:,ii) = theta(:,ii-1) + k*dtheta;
    thetad = theta(:,ii);
    sig(:,ii) = thetad(1:sN); E_GN(:,ii) = thetad(sN+1:end);
    
    %%% projection negative nonphysical parameters %%%
    negind = find(sig(:,ii) < minval);
    sig(negind,ii) = minval;
    
    negind = find(E_GN(:,ii) < minval);
    E_GN(negind,ii) = minval;
    
    %%% resimulate ERT data %%%
    imgGN = mk_image(fmdlC, sig(:,ii));
    Uref = fwd_solve(imgGN);
    Urefel = Uref.meas;
    
    %%% resimulate QSEI data %%%
    [usim,~] = FMDL(E_GN(:,ii),nu,th,constraint,TriQSEIinv,gQSEIinv,cond,force);
    
    %%% recompute cost fun %%%
    Fnew(ii)= norm(um - usim)^2 + norm(Uel - Urefel)^2 + alpha*sum(sqrt((L_TV*(alpha_E*E_GN(:,ii))).^2 + (L_TV*sig(:,ii)).^2)) + lambda2*norm(E_GN(:,ii),1) + lambda1*norm(sig(:,ii),1);
       
    %%%%%%%%%%%%%%%%
    %%% PLOT ERT %%%
    %%%%%%%%%%%%%%%%
    figure(6),
    F=scatteredInterpolant(ExQSEIinv,EyQSEIinv,E_GN(:,ii),'linear');
    int=F(gQSEIinv(:,1),gQSEIinv(:,2));
    trisurf(TriQSEIinv(IO,:),TriQSEIinv.Points(:,1), TriQSEIinv.Points(:,2),int),view(2),colormap('hot'),set(0,'defaulttextInterpreter','latex'),daspect([1 1 1]),
    colorbar('eastoutside');axis tight, title('$E_{inv}$','FontSize',14), axis off, box on;shading interp;
    
    %%%%%%%%%%%%%%%%%
    %%% Plot QSEI %%%
    %%%%%%%%%%%%%%%%%
    figure(7),
    img_plot = mk_image(fmdlC,sig(:,ii));
    show_fem(img_plot),colormap('hot');eidors_colourbar(img_plot);title('$E_{inv}$','FontSize',14);
    drawnow
    
    %%% Stopping criteria (if desired, presently stops at 10 iterations)
end



