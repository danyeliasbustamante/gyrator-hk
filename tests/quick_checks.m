% tests/quick_checks.m
% Sanity checks for the SU(2)-based discrete gyrator implementation.
addpath('../src');

tol = 1e-10;
N = 16;

% gamma = 0 -> identity
[HK0, LK0] = gyrator_hk_master(N, 0);
e0 = max(cellfun(@(h,l) norm(h(:)-l(:)), HK0, LK0));
fprintf('gamma=0: max ||HK - LK||_2 = %.3e (<= %.1e expected)\n', e0, tol);

% gamma = pi/2 -> (nx,ny) <-> (ny,nx) up to a known phase
% Here we just check energy preservation and shape equality of magnitudes.
[~, LKpi2] = gyrator_hk_master(N, pi/2);
mag_ok = all(cellfun(@(l) abs(norm(l(:)) - 1) <= 1e-12, LKpi2));
fprintf('gamma=pi/2: unit-norm check = %d\n', mag_ok);

% gamma = pi -> sign and symmetries (magnitude preserved)
[~, LKpi] = gyrator_hk_master(N, pi);
mag_ok2 = all(cellfun(@(l) abs(norm(l(:)) - 1) <= 1e-12, LKpi));
fprintf('gamma=pi: unit-norm check = %d\n', mag_ok2);

if e0 <= tol && mag_ok && mag_ok2
  fprintf('All quick checks passed.\n');
else
  fprintf('Some checks failed. Review implementation or tolerances.\n');
end
