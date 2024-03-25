function [s_d,s_n] = testing_compute_frictionbem(rcv,bcindex,bc_left,bc_right,f)

% density to calculate increase in pressure with depth
rhog = 2;

% compute stress and displacement kernels for BVP
obs = rcv.xc;

tic
[Kdd,Kdn,Knd,Knn] = geometry.computeTractionKernels(rcv,rcv);

[Gdx,Gdz,Gnx,Gnz] = geometry.computeDisplacementKernels(rcv,obs);

Islip = eye(rcv.N);

toc

bc_dirichlet = find(bcindex==0);

BC_1 = zeros(rcv.N,1);% mixed boundary conditions for 1-component (d/x)
BC_2 = zeros(rcv.N,1);% mixed boundary conditions for 2-component (n/z)

Kd1 = Kdd;
Kd2 = Kdn;
Kn1 = Knd;
Kn2 = Knn;

% all boundary elements that have displacement BCs
ind = false(rcv.N,1);
for i = 1:length(bc_dirichlet)
    ind = ind|rcv.Vpl==bc_dirichlet(i);
end

Kd1(ind,:) = Gdx(ind,:);
Kd2(ind,:) = Gdz(ind,:);
Kn1(ind,:) = Gnx(ind,:);
Kn2(ind,:) = Gnz(ind,:);

% interior (add friction kernels)
ind = rcv.Vpl == 4;
Kd1(ind,:) = Kdd(ind,:)-f*Kdn(ind,:);
Kn1(ind,:) = Knd(ind,:)-f*Knn(ind,:);
Kd2(ind,:) = Islip(ind,:).*0;
Kn2(ind,:) = Islip(ind,:);

%% apply boundary conditions
% left
ind = rcv.Vpl==1;
BC_1(ind) = bc_left(1);
BC_2(ind) = bc_left(2);

% bottom
ind = rcv.Vpl==2;
BC_1(ind) = 0;
BC_2(ind) = 0;

% right
ind = rcv.Vpl==3;
BC_1(ind) = bc_right(1);
BC_2(ind) = bc_right(2);

% interior
ind = rcv.Vpl==4;
BC_1(ind) = 0;%f*rhog*rcv.xc(ind,2);
BC_2(ind) = 0;

BCvec = [BC_1;BC_2];

% construct giant stress kernel matrix
bigK = [Kd1,Kn1;...
        Kd2,Kn2];

%% compute bem solution
% solution vector - [S_d;S_n];
% friction condition on interior fault
% interior
ind = rcv.Vpl == 4;

% tau_d = [Kdd(ind,:),Knd(ind,:)]*[s_d;s_n]
% tau_n = [Kdn(ind,:),Knn(ind,:)]*[s_d;s_n]
% K_tau_d = [Kdd(ind,:),Knd(ind,:)];
% K_tau_n = [Kdn(ind,:),Knn(ind,:)];

% we want to impose:
% |tau_d| <= f|tau_n|
% let us solve the BEM problem naively
disp(['Condition number: ' num2str(cond(bigK),'%.2d')])
solvec = pinv(bigK)*BCvec;
s_d = solvec(1:rcv.N);
s_n = solvec(rcv.N+1:2*rcv.N);

% use naive solution to iteratively add the friction BC
% |K_tau_d*solvec + Kdd.δs_d| <= f|K_tau_n*solvec + ρgz + Kdn.δs_d|
% |tau_d_0 + Kdd.δs_d| <= f|tau_n_0 + Kdn.δs_d|
% let's start the solve by neglecting (Kdn.δs_d)
% if tau_d_0 > 0,
% the KKT conditions are:
% (1) Kdd.δs_d <= f|tau_n_0| - tau_d_0
% (2) tau_d_0 + Kdd.δs_d >= 0 
% (3) (δs_d)(Kdd.δs_d - f|tau_n_0| + tau_d_0) = 0 (complementarity condition)

K = Kdd(ind,ind);
tau_d = [Kdd(ind,:),Knd(ind,:)]*[s_d;s_n];
tau_n = [Kdn(ind,:),Knn(ind,:)]*[s_d;s_n];
p = 3*rcv.xc(ind,2);
sigma_n = abs(p + tau_n);

index = tau_d>=0;
% (tau_d(index) + K(index,:).δs_d) <= f*sigma_n
% (1) tau_d(index) + K(index,:).δs_d - f*sigma_n <= 0
% rewrite as -K(index,:).z + (f*sigma_n) - tau_d(index) >= 0
% 
% for tau_d<0
% -(tau_d(~index) + K(~index,:).δs_d) <= f*sigma_n
% (2) (tau_d(~index) + K(~index,:).δs_d) + f*sigma_n >= 0
% rewrite as K(~index,:).z + (f*sigma_n + tau_d(~index)) >= 0
% 
% problem is how to work with opposite signs of z?
q = [f*sigma_n(index)-tau_d(index);...
     f*sigma_n(~index)+tau_d(~index)];
M = [-K(index,:);K(~index,:)];
% q = f*sigma_n(~index)+tau_d(~index);
% M = K(~index,~index);

[w,z,retcode] = LCPSolve(M,q);


end