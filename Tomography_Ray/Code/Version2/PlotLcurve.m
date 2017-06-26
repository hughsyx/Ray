function PlotLcurve( sm, Idx_T, lambda, rrho, deta, idx_c )
% Select and plot model based on relative residual error tolerance
% Input:
% sm: singular value of G matrix
% Idx_T: the index to truncate G's singular values
% lambda: trade-off parameter in inversion
% rrho: relative residual error vector (norm(d-G*m)/norm(d))
% deta: model change norm vector (norm(dm))
% idx_c: index for selected lambda

        lambda_c = lambda(idx_c);
        rrho_c = rrho(idx_c);
        deta_c = deta(idx_c);
        figure; 
        subplot(222)
        loglog(rrho,deta,'b.-');
        hold on; loglog(rrho_c,deta_c,'r+','markersize',10);axis tight;
        title({['L-curve of ',num2str(length(lambda)),' different \lambda values '],...
        ['using LSQR to minimize || G\deltam - \deltad ||^2_2 + \lambda^2 || \deltam ||^2_2 '],...
        ['( \deltam = m - m0, \deltad = d - d0, lambda chosen = ',num2str(lambda_c),' ) ']},'fontsize',16);
        xlabel('|| ( d - G*m ) ||_2 / || d ||_2');
        ylabel('|| \deltam ||_2');
        axis([min(rrho) max(rrho) min(deta) max(deta)]);
        subplot(221)
        loglog(lambda,deta,'b.-');
        hold on; loglog(lambda_c,deta_c,'r+','markersize',10);
        set(gca,'XDir','reverse');
        title({'Curve of || \deltam ||_2 with different lambda ',...
            ['( for chosen lambda, || \deltam ||_2 = ',num2str(deta_c),' ) ']});
        xlabel(' \lambda ');
        ylabel('|| \deltam ||_2');
        axis([min(lambda) max(lambda) min(deta) max(deta)]);
        subplot(224)
        loglog(rrho,lambda,'b.-');
        hold on; loglog(rrho_c,lambda_c,'r+','markersize',10);
        set(gca,'YDir','reverse');
        title({'Curve of || ( d - G*m ) ||_2 / || d ||_2 with different lambda ',...
            ['( for chosen lambda, || ( d - G*m ) ||_2 / || d ||_2 = ',num2str(rrho_c),' ) ']});
        ylabel(' \lambda ');
        xlabel('|| ( d - G*m ) ||_2 / || d ||_2');
        axis([min(rrho) max(rrho) min(lambda) max(lambda)]);
        subplot(223)
        sm_c=sm(Idx_T);
        loglog(sm,'b.-');hold on;loglog(Idx_T,sm_c,'g+','markersize',10);
        axis tight;
        xlabel('No. of singular values of G ');
        ylabel('Singular values ');
        title({'Singular values "s" of G ',['(cond(G) = ',num2str(max(sm)/min(sm)),...
            ' max(s) = ',num2str(max(sm)),' , min(s) = ',num2str(min(sm)),' ) '],...
            ['( truncate singular values at No. ',num2str(Idx_T),' where s = ',num2str(sm_c),' )']});
        
        
end

