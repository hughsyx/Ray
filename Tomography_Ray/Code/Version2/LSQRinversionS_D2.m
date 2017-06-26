%%% This is to do tomographic inversion with regularized LSQR (Simple Version) %%%
%%% Shaoyu Lu July 12th @ WesternGeco                                          %%%

function [m_c,lambda_c,InvTime1,RMS]=LSQRinversionS_D2(plotLcurve,G,d,k0,k1,xi,zi,Nsrc,Nrec,frzx,lambda_s)
%%%% D for using Differencial Matrix to smooth result
%%%% frzx is the ratio of zflat/xflat
if (nargin==7)
    Nsrc=1;Nrec=length(d);frzx=1;lambda_s=-10000;
end
if (nargin==9)
    frzx=1;lambda_s=-10000;
end
if (nargin==10)
    lambda_s=-10000;
end
nxi=length(xi);
nzi=length(zi);
Dx1=diag([ones(nxi-2,1);0;0])-2*diag([ones(nxi-2,1);0],1)+diag(ones(nxi-2,1),2);
Dx1=Dx1(1:nxi-2,:);
Dx=kron(Dx1,eye(nzi));
Dz1=diag([ones(nzi-2,1);0;0])-2*diag([ones(nzi-2,1);0],1)+diag(ones(nzi-2,1),2);
Dz1=Dz1(1:nzi-2,:);
Dz=kron(eye(nxi),Dz1);
D=[1/(1+frzx)*Dx;frzx/(1+frzx)*Dz];
%% solve with lsqr_b for initial model
if k0<1
    m0=zeros(size(G,2),1);
else
% perform LSQR iteration without reorthogonalization, and plot the
% corresponding discrete L-curve
[MM]=lsqr_b(G,d,k0);
% This "lsqr_b.m" performs k steps of the LSQR Lanczos bidiagonalization 
% algorithm applied to the system    min || A x - b || .
% The routine returns all k solutions, stored as columns of
% the matrix X.  The solution norm and residual norm are returned
% in eta and rho, respectively.

% figure; 
% plot_lc(rho_lsqr_b,eta_lsqr_b,'o');
% title({['L-curve of ',num2str(k_lsqr_b),' iterations '],...
%     ['using LSQR Lanczos bidiagonalization algorithm ']});

m0=MM(:,k0);
Niter_lsqr_b=size(MM,2);rho0_iter=zeros(Niter_lsqr_b,1);
for iter=1:Niter_lsqr_b
    rho0_iter(iter)=norm(d-G*MM(:,iter));
end
end

%% solve with splsqr for updated model

% perform subspace "preconditioned" LSQR iteration without reorthogonalization 
% including computation of the filter factors, and display the filter
   
[U,s,V] = csvd(G); % find singular value of A
sm=s;
%%%% set a minimum singular value to truncake singular value
sm_tol=s(1)/1000; % set smallest singular value tolerance
[junk,IXs]=sort(abs(sm-sm_tol),'ascend');
Idx_T=IXs(1);sm_c=sm(Idx_T);
% fprintf(['\nminimum singular value of G is ',num2str(sm(end)),'\n\n']);
% figure;loglog(sm);
% hold on;loglog(Idx_T,sm_c,'r+','markersize',10);
% xlabel('No. of singular values of G ');
% ylabel('Singular values ');
% title({'Singular values "s" of G ',['(cond(G) = ',num2str(max(sm)/min(sm)),...
%     ' max(s) = ',num2str(max(sm)),' , min(s) = ',num2str(min(sm)),' ) '],...
%     ['( truncate singular values at No. ',num2str(Idx_T),' where s = ',num2str(sm_c),' )']});
sm1=zeros(size(sm));
sm1(1:Idx_T)=sm(1:Idx_T);
sm1(Idx_T+1:end)=ones(length(sm1)-Idx_T,1)*sm(Idx_T);
sm=sm1;

%%% produce L-curve to find appropriate regularization parameter lambda
if lambda_s>=0
    lambda_lc=lambda_s;
    tic
    [xl]=lsqr_b([G;lambda_lc*D],[d-G*m0;zeros(size(D,1),1)],k1);
    %[xl]=splsqr(G,d-G*m0,lambda_lc,eye(size(G,2))); % use default maxit
    toc
    InvTime1=toc;
    fprintf(['\nTime costs for one inversion is ',num2str(InvTime1),' seconds\n']);
    Niter=size(xl,2); rhoc_iter=zeros(Niter,1);
    for iter=1:Niter
        rhoc_iter(iter)=norm(d-G*m0-G*xl(:,iter));
    end
    dml=xl(:,end);
else
    [ dml,lambda_lc,rho_l,eta_l,lambda_l,rhoc_iter,InvTime1 ] = SelectModel_Lcurve( G, U, sm, d-G*m0, k1, plotLcurve );
        % find the index of lambda_lc in lambda_l
    [junk,IX]=sort(abs(lambda_l-lambda_lc),'ascend');
    idx_l=IX(1);
    %PlotLcurve( s, Idx_T, lambda_l, rho_l, eta_l, idx_l );
end

% for ils=1:length(lambda_la)
%     [junk,IX]=sort(abs(lambda_l-lambda_la(ils)),'ascend');
%     idx_la(ils)=IX(1);
%     m_a(:,ils)=dma(:,ils)+m0;
%     PlotInversion( m_a(:,ils),lambda_la(ils),m0,xi,zi,G,d,rho0_iter,rhoc_iter,Nsrc,Nrec )
% end
%%% compute final results
dm_c=dml; lambda_c=lambda_lc;
m_c=dm_c+m0;
%%% plot results
%PlotInversion( m_c,lambda_c,m0,xi,zi,G,d,rho0_iter,rhoc_iter,Nsrc,Nrec )
RMS = rhoc_iter(end);
end