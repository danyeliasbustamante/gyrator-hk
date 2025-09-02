function [HK, LK, SHW, mgrid] = gyrator_hk_master(N, gamma)
% GYRATOR_HK_MASTER  Unifies and runs the n=2 (pi/4) demo in a single file.
%   [HK, LK, SHW, mgrid] = gyrator_hk_master(N, gamma)
%   - N: lattice size (even). Default N = 16 (j = 8).
%   - gamma: gyration angle. Default gamma = pi/4.
%
% Returns:
%   HK   : cell{3} with n=2 base modes (real and orthonormal)
%   LK   : cell{3} with rotated modes (including the paper's global phase)
%   SHW  : cell{3} with the part to plot (Re if μ>=0, Im if μ<0)
%   mgrid: struct with mx,my (grid from -j:j)
%
% This file contains all dependencies as local subfunctions:
%   - gyrator_mode_eq27
%   - hermiteKravchuk2D
%   - kravchukOrthonormal
%   - wignerLittleD
%
% Quick usage in MATLAB:
%   gyrator_hk_master;                   % run with N=16, gamma=pi/4
%   [HK,LK,SHW,~] = gyrator_hk_master(32, pi/6);  % larger grid and another angle
%
% -------------------------------------------------------------------------

    if nargin < 1 || isempty(N),     N = 16;      end
    if nargin < 2 || isempty(gamma), gamma = pi/4; end
    assert(mod(N,2)==0 && N>0, 'N must be even and positive.');

    % --- Grid ---
    j = N/2;
    [mx,my] = meshgrid(-j:j);
    mgrid.mx = mx;  mgrid.my = my;

    % --- Pairs (nx,ny) with n=2: μ = +2, 0, -2 ---
    pares  = [2 0; 1 1; 0 2];
    K = size(pares,1);

    % --- n=2 HK basis (to align phase when μ=0) ---
    B20 = hermiteKravchuk2D(2,0,mgrid);
    B11 = hermiteKravchuk2D(1,1,mgrid);
    B02 = hermiteKravchuk2D(0,2,mgrid);
    Bmat = [B20(:), B11(:), B02(:)];   % orthonormal columns

    % Containers
    HK  = cell(K,1);      % original HK modes (real)
    LK  = cell(K,1);      % rotated modes (with the paper's global phase)
    SHW = cell(K,1);      % what we plot (Re if μ>=0, Im if μ<0)

    % ---------- Compute the 3 modes ----------
    for k = 1:K
        nx = pares(k,1); ny = pares(k,2);
        mu = nx - ny;

        % Base HK mode
        HK{k} = hermiteKravchuk2D(nx, ny, mgrid);

        % Gyration and global phase (paper): Λ = e^{+iπμ/2} · G_{π/4}{φ}
        phi_rot = gyrator_mode_eq27(nx, ny, gamma, mgrid);
        LK{k}   = exp(1i*pi*mu/2) * phi_rot;

        % --- Phase alignment for μ = 0 (make the mode real) ---
        if mu == 0
            c     = Bmat' * LK{k}(:);          % [c20; c11; c02] in the HK basis
            theta = angle(c(1) + c(3));        % enforce c20,c02 to be real
            LK{k} = LK{k} * exp(-1i*theta);    % apply the proper global phase
        end

        % Which part to display
        if mu >= 0
            SHW{k} = real(LK{k});
        else
            SHW{k} = imag(LK{k});
        end
    end

    % ---------- 3×3 figure ----------
    figure('Color','w');
    tiledlayout(3,3,'TileSpacing','compact','Padding','compact')

    for k = 1:K
        nx = pares(k,1); ny = pares(k,2); mu = nx - ny;

        % Column 1: HK
        nexttile
        imagesc(HK{k}); axis square; set(gca,'YDir','normal'); colorbar
        c = max(abs(HK{k}(:))); if c==0, c=1; end, caxis([-c c])
        title(sprintf('\\phi_{%d,%d} (HK)',nx,ny),'Interpreter','tex')

        % Column 2: |Λ|
        nexttile
        imagesc(abs(LK{k})); axis square; set(gca,'YDir','normal'); colorbar
        title('|\\Lambda_{n=2,\\mu}|, \\gamma=\\pi/4','Interpreter','tex')

        % Column 3: Re/Im depending on the sign of μ
        nexttile
        imagesc(SHW{k}); axis square; set(gca,'YDir','normal'); colorbar
        c = max(abs(SHW{k}(:))); if c==0, c=1; end, caxis([-c c])
        if mu >= 0
            ttl = sprintf('Re{\\Lambda} (\\mu=%+d)', mu);
        else
            ttl = sprintf('Im{\\Lambda} (\\mu=%+d)', mu);
        end
        title(ttl,'Interpreter','tex')
    end

    colormap gray
    sgtitle('n=2 \\rightarrow HK modes and rotation to \\gamma=\\pi/4 (\\mu=0)', ...
            'Interpreter','tex')

end % === end of main function ===


% ======================================================================
%  Local subfunctions (kept in this file)
% ======================================================================

function phi_rot = gyrator_mode_eq27(nx, ny, gamma, mgrid)
% φ_{nx,ny}^{(γ)}(mx,my) = ∑_{nx'+ny'=n} d^{j_spin}_{μ/2, μ'/2}(2γ) ·
%                          e^{-iπμ/4} e^{+iπμ'/4} · φ_{nx',ny'}(mx,my)
%
% Requires:  wignerLittleD.m   (Wigner small-d coefficient)
%            hermiteKravchuk2D.m  (2D base modes, orthonormal)
%
% Phase notes (convention used here):
%   For γ = π/2 (β = π) we have
%       φ^{(π/2)}_{nx,ny} = (-1)^{nx} e^{-iπ (nx - ny)/2} · φ_{ny,nx}
%   which is equivalent to
%       φ^{(π/2)}_{nx,ny} = (-1)^{ny} e^{+iπ (nx - ny)/2} · φ_{ny,nx}
%   since μ = nx - ny and (-1)^{nx} e^{-iπμ/2} = (-1)^{ny} e^{+iπμ/2}.

    % --------- Validation / parameters ----------
    N = size(mgrid.mx,1) - 1;              % N = 2*j_grid
    assert(all(size(mgrid.mx)==size(mgrid.my)) && size(mgrid.mx,2)==N+1, ...
        'mgrid.mx/my must be (N+1)x(N+1) with N = size(mx,1)-1.');
    assert(nx>=0 && ny>=0 && mod(nx,1)==0 && mod(ny,1)==0, ...
        'nx, ny must be nonnegative integers.');
    assert(nx+ny <= N, 'nx+ny must be <= N.');

    % --------- Invariant subspace (fixed n) ------
    n       = nx + ny;          % total energy
    j_spin  = n / 2;            % SU(2) spin
    mu      = nx - ny;          % input magnetic number

    % --------- Accumulator -----------------------
    phi_rot = zeros(size(mgrid.mx));   % becomes complex due to phases

    % --------- Sum over (nx',ny') with n'=n ----
    for nxp = 0:n
        nyp = n - nxp;
        mup = nxp - nyp;

        % small-d coefficient (β = 2γ)
        dcoeff = wignerLittleD(j_spin, mu/2, mup/2, 2*gamma);

        % geometric phase (convention: e^{-iπμ/4} e^{+iπμ'/4})
        phase  = exp(-1i*pi*mu/4) * exp( 1i*pi*mup/4 );

        % base mode
        phi_term = hermiteKravchuk2D(nxp, nyp, mgrid);

        % accumulation
        phi_rot = phi_rot + phase * dcoeff * phi_term;
    end

    % --------- (Optional) Renormalize ----------
    % The gyration is unitary; rounding can deviate ~1e-15.
    if abs(1 - norm(phi_rot(:))) > 1e-12
        phi_rot = phi_rot / norm(phi_rot(:));
    end
end


function phi = hermiteKravchuk2D(nx, ny, mgrid)
    % 1) parameters
    N      = size(mgrid.mx,1)-1;  
    j_grid = N/2;
    mVec   = -j_grid:j_grid;

    % 2) orthonormal polynomials
    Kx = kravchukOrthonormal(nx, mVec, N);
    Ky = kravchukOrthonormal(ny, mVec, N);

    % 3) MATLAB indices
    idxX = mgrid.mx + j_grid + 1;
    idxY = mgrid.my + j_grid + 1;

    % 4) Mode WITHOUT weight
    phi0 = Ky(idxY) .* Kx(idxX);

    % 5) Euclidean normalization
    phi  = phi0 / norm(phi0(:));
end


function K = kravchukOrthonormal(n, mVec, N)
    % Kravchuk polynomials orthonormal in Euclidean norm (unit norm)
    j_grid = N/2;
    x      = (mVec + j_grid)';        % column 0,1,...,N

    % binomial weight p=1/2:
    w = gamma(N+1) ./ ( gamma(x+1).*gamma(N-x+1) ) * (1/2)^N;

    % build Kraw:
    Kraw = zeros(size(x));
    for j0 = 0:n
        c1 = gamma(x+1)    ./ ( gamma(j0+1)    .* gamma(x - j0 +1) );
        c2 = gamma(N-x+1)  ./ ( gamma(n-j0+1)  .* gamma((N-x) - (n-j0) +1) );
        Kraw = Kraw + (-1)^j0 * (c1 .* c2);
    end

    % apply sqrt(w) and normalize:
    K = (Kraw .* sqrt(w))';
    K = K / norm(K);
end


function d = wignerLittleD(j, m, mp, beta)
    % d^j_{m,mp}(beta), j integer or half-integer; m, mp with same parity
    tol = 1e-12;
    b = mod(beta, 2*pi);

    % Stable limiting cases
    if abs(b) < tol
        d = double(abs(m-mp) < tol); % Kronecker δ_{m,mp}
        return
    elseif abs(b - pi) < tol
        % d^j_{m,mp}(π) = (-1)^{j-m} δ_{m,-mp}
        d = ((-1)^(j - m)) * double(abs(mp + m) < tol);
        return
    end

    G = @(x) gamma(x+1);
    cb = cos(beta/2);  sb = sin(beta/2);

    % Prefactor and correct summation limits
    pref = sqrt( G(j+m) * G(j-m) * G(j+mp) * G(j-mp) );
    kmin = max(0, m - mp);
    kmax = min(j + m, j - mp);

    dsum = 0;
    for k = kmin:kmax
        num = (-1)^(k - mp + m) * cb^(2*j + m - mp - 2*k) * sb^(mp - m + 2*k);
        den = G(j + m - k) * G(k) * G(j - mp - k) * G(k - m + mp);
        dsum = dsum + pref * num / den;
    end
    d = dsum;
end
