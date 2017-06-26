function PlotInversion( m_c,lambda_c,m0,xi,zi,G,d,rho0_iter,rhoc_iter,Nsrc,Nrec )
%Plot inversion result
%   obtain initial model by min{||Ax-b||}, update by
%   min{||A*dx-db||^2+lambda^2*||dx||^2}

nxi=length(xi);nzi=length(zi);
dd0=d-G*m0; dd_c=d-G*m_c;
Niter0=length(rho0_iter); Niter_c=length(rhoc_iter);
mscale=[min(min(m_c(:)),min(m0(:))),max(max(m_c(:)),max(m0(:)))];
ddscale=[min(min(dd0(:)),min(dd_c(:))),max(max(dd0(:)),max(dd_c(:)))];

figure;
subplot(231);
plot([1:Niter0+Niter_c]',[rho0_iter;rhoc_iter],'b.-');
hold on; plot(Niter0,rho0_iter(end),'r+','markersize',10);axis tight;
hold on; text(Niter0,rho0_iter(end),'\leftarrow m_0 ');
xlabel('Iteration times ');ylabel('Data residual norm ');
title({'Data residual changes with iteration times ',...
    ['Final data residual = ',num2str(rhoc_iter(end)),' ']});
subplot(232);
imagesc(xi,zi,reshape(m0,nzi,nxi),mscale);axis image; colorbar;
xlabel('x ');ylabel('z ');
title('Initial model obtained by min(||d-G*m||) ');
subplot(233);
imagesc(xi,zi,reshape(m_c,nzi,nxi),mscale);axis image; colorbar;
xlabel('x ');ylabel('z ');
title({'Final model updated by min(||\deltad-G*\deltam||^2+\lambda*||\deltam||^2) ',...
    ['( \lambda = ',num2str(lambda_c),' ) ']});
if Nsrc==1 %one dimension data
    subplot(234);
    plot(d,'r.-');
    xlabel('No. of data '); ylabel('Data ');
    title('Observed data values ');
    subplot(235);
    plot(dd0,'b.-');
    xlabel('No. of data '); ylabel('Data residuals ');
    title('Data residuals computed using initial model ');
    subplot(236);
    plot(dd_c,'b.-');
    xlabel('No. of data '); ylabel('Data residuals ');
    title('Data residuals computed using updated model ');
else
    subplot(234);
    imagesc(reshape(d,Nrec,Nsrc));
    axis image; colorbar;
    xlabel('No. of sources '); ylabel('No. of receivers ');
    title('Observed data values ');
    subplot(235);
    imagesc(reshape(dd0,Nrec,Nsrc),ddscale);
    axis image; colorbar;
    xlabel('No. of sources '); ylabel('No. of receivers ');
    title('Data residuals computed using initial model ');
    subplot(236);
    imagesc(reshape(dd_c,Nrec,Nsrc),ddscale);
    axis image; colorbar;
    xlabel('No. of sources '); ylabel('No. of receivers ');
    title('Data residuals computed using updated model ');
end
    
end

