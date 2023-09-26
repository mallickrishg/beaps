function [slip_d,slip_n] = compute_bemsolution(rcv,bcindex,bc_left,bc_right)

% compute stress and displacement kernels for BVP
obs = rcv.xc;

tic
[Kdd,Kdn,Knd,Knn] = geometry.computeTractionKernels(rcv,rcv);

[Gdx,Gdz,Gnx,Gnz] = geometry.computeDisplacementKernels(rcv,obs);
toc

bc_dirichlet = find(~bcindex);

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

BCvec = [BC_1;BC_2];

% construct giant stress kernel matrix
bigK = [Kd1,Kn1;...
        Kd2,Kn2];

%% compute bem solution
% solution vector - [S_d;S_n];
solvec = pinv(bigK)*BCvec;
% solvec = regress(BCvec,bigK);
% solvec = lsqlin(bigK,BCvec,[],[],[],[],-10.*ones(rcv.N*2,1),10.*ones(rcv.N*2,1));
% solvec = bigK\BCvec;

slip_d = solvec(1:rcv.N);
slip_n = solvec(rcv.N+1:2*rcv.N);

end