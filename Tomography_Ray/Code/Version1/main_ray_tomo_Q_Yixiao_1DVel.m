%invert Q
%% set the parameters
alpha_init = 0.00005*ones(nz,nx); % initial guess of alpha
ngrid = size(G_matrix,2); ntrace = size(G_matrix,1);
lambda2 = lambda/10;  % no idea if this is a proper one
%% input observed frequency shift freq_shift_obs(ntrace,1), sigma_f, Ind_use
%sigma_f = 10.3154; %standard deviation of input signal
%freq_shift_obs = rand(ntrace,1);

load([Tomo_path, '/Models/cross_well_att_f_obs.mat']);
sigma_f = sqrt(219);%f_sigma_initial*1.5;
freq_shift_obs = fc_shift;
Ind_use = I_use;
clear('f_sigma_initial','fc_shift','I_use');
%% build the new G_matrix and calculate freq shift based on initial guess of alpha
G2_matrix = sparse(ntrace,ngrid+1);
G2_matrix(:,1:end-1) = G_matrix;
%G2_matrix(:,end) = +1/sigma_f^2*ones(ntrace,1);
G2_matrix(:,end) = 1/ones(ntrace,1);
freq_shift_guess = G_matrix*alpha_init(:);

%% inversion 
d_freq_shift = (freq_shift_obs)/sigma_f^2 - freq_shift_guess;
        cvx_begin
        cvx_quiet(true)
        variable Alpha_Q(ngrid+1)
        minimize(norm(d_freq_shift(Ind_use) - G2_matrix(Ind_use,:) * Alpha_Q(:)) + lambda * norm(Alpha_Q(1:end-1)))
        subject to
            norm(Alpha_Q(end)) <= 0.1*fc_init
        %minimize(norm(d_freq_shift(Ind_use) - G2_matrix(Ind_use,:) * Alpha_Q(:)) + lambda * norm(Alpha_Q(1:end-1)) + lambda2 * Alpha_Q(end))
        
        %minimize(lambda * norm(Alpha_Q(1:end-1)) + Alpha_Q(end))
        %subject to 
        %    G2_matrix(Ind_use,:) * Alpha_Q(:) == d_freq_shift(Ind_use)
        cvx_end
        
 %delta_f = -Alpha_Q(end);
 delta_f = -Alpha_Q(end)*sigma_f^2;
 dalpha = reshape(Alpha_Q(1:end-1),nz,nx);
 dalpha = modelSmooth(dalpha, 3);
 
 %% get Q
 Q = pi./((dalpha+alpha_init).*V_collec(:,:,end));