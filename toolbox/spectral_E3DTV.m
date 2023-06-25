function img_denoise = spectral_E3DTV(sino, G, x0, par)
%%  Code for paper (reconstruction part):
%
%  D. Li, D. Zeng, S. Li, Y. Ge, Z. Bian, J. Huang, and J. Ma, "MDM-PCCT:
%  Multiple Dynamic Modulations for High-Performance Spectral PCCT Imaging,"
%  in IEEE Transactions on Medical Imaging, vol. 39, no. 11, pp. 3630-3642,
%  Nov. 2020, doi: 10.1109/TMI.2020.3001616.
%
%  ------------------------------------------------------------------------
%
%% Solving the iterative reconstruction model:
%
%             min_{X} = 1/2||A*X-Y||_2^2 + alpha*E3DTV(X)
%                s.t. D_x(X) = U_x*V_x, V_x'*V_x=I
%                     D_y(X) = U_y*V_y, V_y'*V_y=I
%                     D_z(X) = U_z*V_z, V_z'*V_z=I
%                     rank(U_{x,y,z})=r
%     ==============================================================
%        Y is the sinogram data, A is the system matirx, 
%        X is the desired image, alpha is the hypar-parameter,
%        E3DTV(.) is the prior term based on enhanced 3DTV which encodes
%                 the sparse property of the image in the subspace of the
%                 image gradient maps, 
%        D is difference operator
%        U_{x,y,z} are the subspace of the gradient maps,
%        V_{x,y,z} are the transformation matrixes, 
%        T is the transformation operation, 
%        I is the identity matrix, 
%        r is the rank of the multi-energy image.
%
%  ------------------------------------------------------------------------
%
% Input:
%        sino: Input noisy spectral or photon-couting sinogram data [ndetector*nview*nbin]
%        G: System matrix of the imaging system, 2D
%        x0: initialize reconstruction image of the multi-bin, it can be the
%            FBP-reconstructed image
%        par: parameters of the model and the optimization procedure
%
% Output:
%        img_denoise: final reconstructed spectral or PCCT image.
%
%  ------------------------------------------------------------------------
%
% by Danyang Li (email: dyli0730@gmail.com)
% 2019/7/17
% 
%  ------------------------------------------------------------------------
%

%% parameter setting
x    = x0;
y    = sino;
nbin = size(y, 3);

% parameters
alpha      = par.alpha;
v_min      = par.v_min;
v_max      = par.v_max;
recon_iter = par.recon_iter;
tau        = par.tau;
rank       = par.rank;
tol        = 1e-6;
rho        = par.rho;
max_mu     = 1e6;
mu         = par.mu;
lambda     = par.lambda;

%% Initializing optimization variables
% normalize the image
disp(['Initiailzing optimization variables ...'])
x_init = (x - v_min)./(v_max - v_min);

[M,N,p] = size(x_init);
sizeD   = size(x_init);
D       = zeros(M*N,p);
scal    = zeros(p,1);
for i=1:p
    bandp     = x_init(:,:,i);
    scal(i,1) = sqrt(sum(bandp(:).^2)/(M*N));
    D(:,i)    = bandp(:);  % unfolding on the spectral mode
end
normD   = norm(D,'fro');

R              = randn(M*N,p);
E              = zeros(M*N,p);
% U_x and V_x initial
tv_x           = diff_x(R,sizeD);
tv_x           = reshape(tv_x,[M*N,p]);
[U_x,S_x,V_x]  = svd(tv_x,'econ');
U_x            = U_x(:,1:rank(1))*S_x(1:rank(1),1:rank(1));
V_x            = V_x(:,1:rank(1));
% U_y and V_y initial
tv_y           = diff_y(R,sizeD);
tv_y           = reshape(tv_y,[M*N,p]);
[U_y,S_y,V_y]  = svd(tv_y,'econ');
U_y            = U_y(:,1:rank(2))*S_y(1:rank(2),1:rank(2));
V_y            = V_y(:,1:rank(2));
% tv_z initial
tv_z           = diff_z(R,sizeD);
tv_z           = reshape(tv_z,[M*N,p]);
[U_z,S_z,V_z]  = svd(tv_z,'econ');
U_z            = U_z(:,1:rank(3))*S_z(1:rank(3),1:rank(3));
V_z            = V_z(:,1:rank(3));
Mul1 =zeros(size(D));  % multiplier for D-X-E
Mul2 =zeros(size(D));  % multiplier for Dx_X-U_x*V_x
Mul3 =zeros(size(D));  % multiplier for Dy_X-U_y*V_y
Mul4 =zeros(size(D));  % multiplier for Dz_X-U_z*V_z

%% main loop
for iter = 1:recon_iter
    disp(['Iteration: ', num2str(iter), '/', num2str(recon_iter), ' ...'])

    %% -Updata R (ZZ)
    disp(['Update Z in iteration ', num2str(iter), 'th ...'])
    r_t           = R(:);
    Cha           = D-E;
    Mu_x          = U_x*V_x';
    Mu_y          = U_y*V_y';
    Mu_z          = U_z*V_z';
    r_t           = myPCG_sstv(r_t,Cha,Mu_x,Mu_y,Mu_z,Mul1,Mul2,Mul3,Mul4,mu,sizeD);
    R             = reshape(r_t,[M*N,p]);
    for i = 1:p
        post_scal = sqrt(sum(R(:,i).^2)/(M*N));
        R(:,i) = R(:,i) * scal(i) ./ post_scal; % scal correction
    end

    ZZ = reshape(R,[M,N,p]);
    ZZ(ZZ<0) = 0;
    % inv normalize
    ZZ = ZZ.*(v_max - v_min) + v_min;

    %% -Update U_x and U_y and U_z
    disp(['Update U in iteration ', num2str(iter), 'th ...'])
    tmp_x         = reshape(diff_x(R,sizeD),[M*N,p]);
    tmp_x         = tmp_x+Mul2/mu;
    U_x           = softthre(tmp_x*V_x, tau/mu);
    tmp_y         = reshape(diff_y(R,sizeD),[M*N,p]);
    tmp_y         = tmp_y+Mul3/mu;
    U_y           = softthre(tmp_y*V_y, tau/mu);
    tmp_z         = reshape(diff_z(R,sizeD),[M*N,p]);
    tmp_z         = tmp_z+Mul4/mu;
    U_z           = softthre(tmp_z*V_z, tau/mu);

    %% -Update V_x and V_y and V_z
    disp(['Update V in iteration ', num2str(iter), 'th ...'])
    [u,~,v]       = svd(tmp_x'*U_x,'econ');
    V_x           = u*v';
    [u,~,v]       = svd(tmp_y'*U_y,'econ');
    V_y           = u*v';
    [u,~,v]       = svd(tmp_z'*U_z,'econ');
    V_z           = u*v';

    %% -Update E
    disp(['Update E in iteration ', num2str(iter), 'th ...'])
    E             = softthre(D-R+Mul1/mu, lambda/mu);

    %% -Update X in cpu version, it costs too much running time
    disp(['Update X in iteration ', num2str(iter), 'th ...'])
    for ibin = 1:nbin
        x_t     = x(:,:,ibin);
        x_t     = x_t - alpha*(x_t - ZZ(:,:,ibin));
        g       = G'*(y(:,:,ibin) - G*x_t);
        Gg      = G*g;
        lr      = sum(g(:).^2)/sum(Gg(:).^2);
        x(:,:,ibin) = x_t + lr*g;
    end
    
    %% stop criterion
    leq1 = D - R - E;
    leq2 = reshape(diff_x(R,sizeD),[M*N,p])- U_x*V_x';
    leq3 = reshape(diff_y(R,sizeD),[M*N,p])- U_y*V_y';
    leq4 = reshape(diff_z(R,sizeD),[M*N,p])- U_z*V_z';
    stopC1 = norm(leq1,'fro')/normD;
    stopC2 = max(abs(leq2(:)));
    if stopC1<tol && stopC2<tol
        break;
    else
        Mul1 = Mul1 + mu*leq1;
        Mul2 = Mul2 + mu*leq2;
        Mul3 = Mul3 + mu*leq3;
        Mul4 = Mul4 + mu*leq4;
        mu = min(max_mu,mu*rho);
    end
    

end

%% output
img_denoise = x;

end


