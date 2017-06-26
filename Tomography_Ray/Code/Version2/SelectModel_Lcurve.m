function [ dml,lambda_lc,rho,eta,lambda,rhol_Iter,InvTime1 ] = SelectModel_Lcurve( A, U, sm, b, kl, plotLcurve )
% Select and plot model based on relative residual error tolerance
% Input:
% [U,S,V] = csvd(A);
% solving Ax=b by min{||Ax-b||^2+lambda^2*||x||^2)


        [lambda_lc,rho,eta,lambda] = l_curve(plotLcurve,U,sm,b,'Tikh');
        if strcmp(plotLcurve,'yes');
            axis([min(rho) max(rho) min(eta) max(eta)]);
        end
        tic
        [xl]=lsqr_b([A;lambda_lc*eye(size(A,2))],[b;zeros(size(A,2),1)],kl);
        %[xl]=splsqr(A,b,lambda_lc,eye(size(A,2))); % use default maxit
        toc
        InvTime1=toc;
        fprintf(['\nTime costs for one inversion is ',num2str(InvTime1),' seconds\n\n\n']);
        dml=xl(:,end);rhol_Iter=zeros(size(xl,2),1);
        for it=1:size(xl,2)
        rhol_Iter(it)=norm(A*xl(:,it)-b);
        end
        rho=rho/norm(b);
        
        
end

