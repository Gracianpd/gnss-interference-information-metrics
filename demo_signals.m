m = ti_signals();

cfg = struct();
cfg.Fs  = 10000;
cfg.T   = 0.2;
cfg.f0  = 0.02*cfg.Fs;  % 200 Hz
cfg.A   = 1;
cfg.phi = 0;

cfg.D   = 0.25;         % rectangular
cfg.Rc  = cfg.Fs/200;    % Nc=10 muestras/chip
cfg.seed = 1;

S = m.generate_all(cfg);
m.plot_all(S);
