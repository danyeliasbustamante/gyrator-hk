% demos/demo_n2_pi4.m
% Minimal reproducible demo: n=2, gamma=pi/4. Saves a figure in ../figs.
addpath('../src');
[HK, LK, SHW, mgrid] = gyrator_hk_master(16, pi/4); %#ok<ASGLU>
saveas(gcf, '../figs/demo.png');
fprintf('Saved figure to figs/demo.png\n');
