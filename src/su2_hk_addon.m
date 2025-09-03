function su2_hk_addon()
% su2_hk_addon  Tests SU(2) + proyectores de simetría para tu base HK (n=2)
% Requisitos: tener en el path tu archivo gyrator_hk_master.m
%
% Tests:
%  A) Ry(2γ) ≈ gyrator_mode_eq27 (módulo fase global)
%  B) Rx por composición = D(π/2, ψ, −π/2)
%  C) Rz deja fijo μ=0
%  D) Unitariedad y no-conmutatividad
%  E) Proyectores C_n: idempotencia, suma a I y chequeo eigen
%  F) Proyectores D_n (±): idempotencia y suma a I
%
% Nota: los proyectores C_n están formulados de forma segura para n par
% (j entero, como n=2). Para j semientero y nfold impar, usa 2*nfold
% (doble recubrimiento) si lo necesitas en el futuro.

    % ---------- Config ----------
    N     = 16;          % tamaño de retícula (par)
    gamma = pi/6;        % ángulo de gyrator en tu master (β = 2*gamma)
    n     = 2;           % este test usa n=2 (j=1), como en tu demo
    beta  = 2*gamma;

    % ---------- Obtenemos HK/LK del master SIN plot ----------
    [HK, LK] = call_master_silent(N, gamma);

    % Matriz de base B (columnas ortonormales) en el orden μ = +n, 0, −n
    B = [HK{1}(:), HK{2}(:), HK{3}(:)];   % ↔ m = +1,0,−1 para n=2

    % ---------- TEST A: Ry vs gyrator ----------
    Ry  = Ry_hk(n, beta);           % convención de fases compatible con eq.27
    errs = zeros(3,1);
    for k = 1:3
        phi_in  = HK{k};
        phi_su2 = apply_rotation_hk_from_B(phi_in, B, Ry);
        % alineamos fase global con la salida LK{k} de tu master
        theta   = angle( sum(conj(LK{k}(:)).*phi_su2(:)) );
        phi_su2 = phi_su2 * exp(-1i*theta);
        errs(k) = norm(LK{k}(:) - phi_su2(:)) / norm(LK{k}(:));
    end
    fprintf('\n[TEST A] Ry(2γ) vs gyrator: errores relativos por modo = [%.2e %.2e %.2e]\n', errs);
    fprintf('          error máximo = %.2e\n', max(errs));

    % ---------- TEST B: Rx por composición vs Euler ----------
    psi   = 0.37*pi;
    Rx1   = Rx_hk(n, psi);                         % Rz(π/2) Ry(ψ) Rz(−π/2)
    Rx2   = Dj_euler_hk(n, pi/2, psi, -pi/2);      % D(α,β,γ) forma Euler z–y–z
    errRx = norm(Rx1 - Rx2, 'fro');
    fprintf('[TEST B] Rx composición vs Euler: ‖Δ‖_F = %.2e\n', errRx);

    % ---------- TEST C: Rz deja fijo μ=0 ----------
    phi    = 1.234;
    phiRz  = apply_rotation_hk_from_B(HK{2}, B, Rz_hk(n, phi));  % μ=0
    errRz0 = norm(phiRz(:) - HK{2}(:)) / norm(HK{2}(:));
    fprintf('[TEST C] Rz sobre μ=0 (≈0): %.2e\n', errRz0);

    % ---------- TEST D: Unitariedad y no-conmutatividad ----------
    c0 = randn(3,1)+1i*randn(3,1); c0 = c0/norm(c0);
    c1 = Ry*c0; c2 = Rz_hk(n,phi)*c1; c3 = Rx1*c2;
    errUnit = abs(norm(c3) - 1);
    fprintf('[TEST D] Unitariedad (preservación de norma): %.2e\n', errUnit);
    noncomm = norm(Rz_hk(n,phi)*Ry - Ry*Rz_hk(n,phi), 'fro');
    fprintf('[TEST D] No-conmutatividad ‖[Rz,Ry]‖_F: %.2e\n', noncomm);

    % ---------- TEST E: Proyectores C_n ----------
    nfold  = 6;
    I3     = eye(2*(n/2)+1);             % para n=2 → 3x3
    Rstep  = Rz_hk(n, 2*pi/nfold);
    S      = zeros(size(I3));
    maxId  = 0; maxEig = 0;

    for kk = 0:nfold-1
        Pk = projector_Cn_hk(n, nfold, kk);
        maxId  = max(maxId, norm(Pk*Pk - Pk, 'fro'));                              % idempotencia
        maxEig = max(maxEig, norm(Rstep*Pk - exp(1i*2*pi*kk/nfold)*Pk, 'fro'));    % eigen-sectores
        S = S + Pk;
    end
    fprintf('[TEST E] C_n: max ‖P^2−P‖=%.2e,  max eigen-check=%.2e,  ‖ΣP−I‖=%.2e\n', ...
            maxId, maxEig, norm(S - I3, 'fro'));

     % ---------- TEST F: Proyectores D_n (±) corregidos ----------
    Fmir  = mirror_y_hk(n);
    errF2 = norm(Fmir*Fmir - eye(size(Fmir)), 'fro');

    Sdn = zeros(size(I3));   maxIdpm = 0;
    visited = false(1, nfold);
    for kk = 0:nfold-1
        if visited(kk+1), continue; end
        kk_neg = mod(nfold - kk, nfold);
        Pkp = projector_Dn_hk(n, nfold, kk, +1);
        Pkm = projector_Dn_hk(n, nfold, kk, -1);
        maxIdpm = max([maxIdpm, norm(Pkp*Pkp-Pkp,'fro'), norm(Pkm*Pkm-Pkm,'fro')]);
        Sdn = Sdn + Pkp + Pkm;
        visited(kk+1)      = true;
        visited(kk_neg+1)  = true;
    end
    fprintf('[TEST F] D_n (fix): ‖F^2−I‖=%.2e,  max idempotencia ± = %.2e,  ‖ΣP±−I‖=%.2e\n', ...
            errF2, maxIdpm, norm(Sdn - I3, 'fro'));

        % ---------- TEST G: Permutación x<->y ----------
    Pxy   = Pxy_hk(n);
    swap  = [3 2 1];                    % φ20 <-> φ02, φ11 -> φ11
    ePerm = zeros(1,3);
    for k = 1:3
        phi_perm = apply_rotation_hk_from_B(HK{k}, B, Pxy);
        phi_tgt  = HK{swap(k)};
        th       = angle( sum(conj(phi_tgt(:)).*phi_perm(:)) );
        ePerm(k) = norm(phi_perm(:)*exp(-1i*th) - phi_tgt(:))/norm(phi_tgt(:));
    end
    fprintf('[TEST G] x<->y swap (Ry(pi)): max error = %.2e\n', max(ePerm));
    % ---------- DEMO: fracciones de simetría en barrido α ----------
    doPlot  = true;
    nfold   = 6;                          % C6 / D6
    alphas  = linspace(0, pi, 121);       % α en [0,π]
    c0      = B' * HK{1}(:);              % estado de partida (φ_{2,0})
    Fop     = mirror_y_hk(n);

    % precompute proyectores
    Pc  = cell(1,nfold); Pp = cell(1,nfold); Pm = cell(1,nfold);
    for k = 0:nfold-1
        Pc{k+1} = projector_Cn_hk(n, nfold, k);
        Pp{k+1} = projector_Dn_hk(n, nfold, k, +1);
        Pm{k+1} = projector_Dn_hk(n, nfold, k, -1);
    end

    E_C  = zeros(nfold, numel(alphas));
    E_Dp = zeros(nfold, numel(alphas));
    E_Dm = zeros(nfold, numel(alphas));
    Par  = zeros(1, numel(alphas));

    for t = 1:numel(alphas)
        beta = 2*alphas(t);          % SU(2): β = 2α
        c    = Ry_hk(n, beta) * c0;  % coeficientes en base HK
        nn   = norm(c)^2;
        for k = 0:nfold-1
            E_C(k+1,t)  = norm(Pc{k+1} * c)^2  / nn;
            E_Dp(k+1,t) = norm(Pp{k+1} * c)^2 / nn;
            E_Dm(k+1,t) = norm(Pm{k+1} * c)^2 / nn;
        end
        Par(t) = real(c' * Fop * c); % ⟨F⟩ ∈ [-1,1]
    end

    if doPlot
        figure('Color','w');
        subplot(3,1,1); plot(alphas, E_C.');  grid on
        xlabel('\alpha'); ylabel('E_k');      title('C_n fractions (n=6)')
        subplot(3,1,2); plot(alphas, E_Dp.'); grid on
        xlabel('\alpha'); ylabel('E_{k,+}');  title('D_n even (+)')
        subplot(3,1,3); plot(alphas, E_Dm.'); grid on
        xlabel('\alpha'); ylabel('E_{k,-}');  title('D_n odd (−)')

        figure('Color','w');
        plot(alphas, Par, 'LineWidth',1.1); grid on
        xlabel('\alpha'); ylabel('\langle F\rangle'); title('Mirror expectation')
    end

end

% ========================= Helpers de integración =========================
function [HK, LK, SHW, mgrid] = call_master_silent(N, gamma)
% Llama gyrator_hk_master sin que dibuje (compatibilidad con y sin doPlot)
    try
        [HK, LK, SHW, mgrid] = gyrator_hk_master(N, gamma, false);
    catch
        [HK, LK, SHW, mgrid] = gyrator_hk_master(N, gamma);
    end
end

% ====================== Núcleo SU(2) en base HK (n fijo) ==================
function R = Rz_hk(n, phi)
    j = n/2;
    m = (j:-1:-j).';              % orden {+j,...,−j} ↔ columnas B
    R = diag(exp(-1i*phi*m));     % diagonal: solo fases por m
end

function R = Ry_hk(n, beta)
    % Sandwich de fases χ=(π/2)m para compatibilidad con tu eq.27
    j   = n/2;
    m   = (j:-1:-j).';
    d   = littleD_matrix(j, beta, m);   % d^j_{m,m'}(β)
    chi = (pi/2)*m;
    R   = diag(exp(-1i*chi)) * d * diag(exp( 1i*chi));
end

function R = Rx_hk(n, psi)
    % Construcción estándar: Rx = Rz(π/2) Ry(ψ) Rz(−π/2)
    R = Rz_hk(n, pi/2) * Ry_hk(n, psi) * Rz_hk(n, -pi/2);
end

function D = Dj_euler_hk(n, alpha, beta, gamma)
    % D^j(α,β,γ) en la base HK (mismo sandwich de fases)
    j   = n/2;
    m   = (j:-1:-j).';
    d   = littleD_matrix(j, beta, m);
    Dz  = diag(exp(-1i*alpha*m));
    Dg  = diag(exp(-1i*gamma*m));
    chi = (pi/2)*m;
    Dm  = Dz * d * Dg;                       % forma estándar en |j,m>
    D   = diag(exp(-1i*chi)) * Dm * diag(exp( 1i*chi)); % pasar a base HK
end

function phi_out = apply_rotation_hk_from_B(phi_in, B, R)
    % Aplica R (3x3) a un modo φ expresado en la subbase HK (columnas de B)
    c       = B' * phi_in(:);               % coeficientes en base HK
    c_rot   = R * c;                        % rotación en el subespacio
    phi_out = reshape(B * c_rot, size(phi_in));
end

% ========================= Wigner-d (estable y local) =====================
function dM = littleD_matrix(j, beta, m_order)
    L = numel(m_order);
    dM = zeros(L,L);
    for a = 1:L
        for b = 1:L
            dM(a,b) = wignerLittleD_scalar(j, m_order(a), m_order(b), beta);
        end
    end
end

function d = wignerLittleD_scalar(j, m, mp, beta)
    % d^j_{m,mp}(β): casos límite y suma estable
    tol = 1e-12;
    b = mod(beta, 2*pi);

    if abs(b) < tol
        d = double(abs(m - mp) < tol);  return
    elseif abs(b - pi) < tol
        d = ((-1)^(round(j - m))) * double(abs(mp + m) < tol);  return
    end

    G  = @(x) gamma(x+1);
    cb = cos(beta/2);  sb = sin(beta/2);

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

% ======================= Proyectores de simetría ==========================
function Pk = projector_Cn_hk(n, nfold, k)
% Proyector a la clase k de C_n actuando como rotación Rz(2π/n) en la base HK.
% Seguro para n par (j entero). Para j semientero y nfold impar, usar 2*nfold.
    j = n/2;  m = (j:-1:-j).';
    p = 0:nfold-1;
    % U(g^p) = Rz(2π p / nfold) → diagonal en m:
    E = exp(-1i*(2*pi/nfold) * (m * p));           % (#m x nfold)
    chi_conj = exp(-1i*(2*pi/nfold) * (k * p));    % (1 x nfold)
    diagVals = mean( E .* chi_conj, 2 );           % (1/n) Σ χ_k^* U(g^p)
    % Redondeo a {0,1} para robustez numérica:
    Pk = diag(double(abs(diagVals) > 1-1e-12));
end

function F = mirror_y_hk(n)
% Reflexión que envía m → −m con la fase (-1)^(j−m) (real y unitaria)
    j = n/2;  m = (j:-1:-j).';
    L = numel(m); F = zeros(L);
    for a = 1:L
        mm = m(a);
        b  = find(abs(m + mm) < 1e-12, 1);   % índice de −m
        F(a,b) = (-1)^(round(j - mm));       % entero siempre
    end
end

function Pkpm = projector_Dn_hk(n, nfold, k, sgn)
% Proyector correcto para D_n:
% 1) forma el par {k, -k}   2) restringe F a ese par   3) P± = 1/2 (Ppar ± Fpar)
    Pc_k   = projector_Cn_hk(n, nfold, k);
    k_neg  = mod(nfold - k, nfold);
    if k_neg == k
        Ppair = Pc_k;                                    % k=0 (y k=n/2 si aplica)
    else
        Pc_neg = projector_Cn_hk(n, nfold, k_neg);
        Ppair  = Pc_k + Pc_neg;                          % par {k,-k}
    end
    F      = mirror_y_hk(n);
    Fpair  = Ppair * F * Ppair;                          % F restringida al par
    Pkpm   = 0.5 * (Ppair + sgn * Fpair);                % sgn = +1 o -1
    Pkpm   = (Pkpm + Pkpm')/2;                           % limpieza numérica
end

function Pxy = Pxy_hk(n)
    % Swap x<->y en la base HK de nivel n (misma convención de fases HK)
    Pxy = Ry_hk(n, pi);   % equivale a intercambiar nx <-> ny (fase global)
end

