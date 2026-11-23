%% SECCION 2 - DSSS/GNSS: 3 frecuencias, curvas vs JSR
% SNR fijo = 5 dB
% Ruidos: white / pink / brown
% Interferencias: NBI / chirp  (SIN CW)
% + Kurtosis espectral (sobre pmf espectral p_i normalizada)
%
% Métricas:
%   - Delta H_norm (espectral)
%   - KL(p_base || p_interf) (espectral)
%   - MI I(Xr;Y) (tiempo), con Xr = x+w y Y = x+w+j
%   - Kurtosis espectral kappa_p (sobre p_i)
%
% Salida:
%   outputs_seccion2/sec2_DSSS_<itf>_<noise>_SNR_5dB.png
%
% Requiere:
%   ti_signals.m
%   ti_noise.m
%   ti_interference.m

clear; close all; clc;

set(groot,'defaultTextInterpreter','none');
set(groot,'defaultLegendInterpreter','none');
set(groot,'defaultAxesTickLabelInterpreter','none');

%% --- directorio salida robusto (.mlx friendly) ---
baseDir = fileparts(which('ti_signals'));
if isempty(baseDir), baseDir = pwd; end

outDir = fullfile(baseDir, 'outputs_seccion2');
if ~exist(outDir,'dir'), mkdir(outDir); end
fprintf("Guardando figuras en: %s\n", outDir);

%% --- módulos ---
sig = ti_signals();
noi = ti_noise();
itf = ti_interference();

%% --- config global ---
Fs = 10000;
T  = 1.0;

% DSSS config
cfgSig = struct();
cfgSig.Fs   = Fs;
cfgSig.T    = T;
cfgSig.A    = 1;
cfgSig.phi  = 0;
cfgSig.Rc   = Fs/10;     % Nc=10
cfgSig.seed = 1;

% 3 frecuencias (bajo/medio/alto)
f0_list = Fs * [0.005, 0.02, 0.08];   % 50, 200, 800 Hz

% SNR fijo
SNRdB = 5;

% Barrido JSR
JSRdB_list = -30:5:15;

% Ruidos a comparar
noiseTypes = {'white','pink','brown'};

% Interferencias (SIN CW)
itfTypes = {'nbi','chirp'};

%% --- interferencias: params fijos salvo JSR/seed ---
cfgItf = struct();
cfgItf.Fs     = Fs;
cfgItf.T      = T;
cfgItf.useDSP = true;
cfgItf.phi    = 0;

% NBI
cfgItf.fj = 2000;
cfgItf.Bj = 400;

% Chirp
cfgItf.f1 = 500;
cfgItf.f2 = 4500;

%% --- Welch (pmf espectral) ---
welch = struct();
welch.winLen  = 4096;
welch.overlap = 2048;
welch.nfft    = 4096;

%% --- MI hist2 ---
miCfg = struct();
miCfg.B     = 32;
miCfg.alpha = 1e-6;
miCfg.L     = 6;   % rango [-L,L] después de estandarizar

%% ==========================================================
% Para cada interferencia y cada tipo de ruido:
%   figura con 4 paneles: DeltaH, KL, MI, Kurtosis
%   En cada panel: 3 curvas (una por f0)
%% ==========================================================

for it = 1:numel(itfTypes)
    itfType = itfTypes{it};

    for nn = 1:numel(noiseTypes)
        noiseType = noiseTypes{nn};

        fig = figure('Name', sprintf('Sec2 DSSS | itf=%s | Ruido=%s | SNR=%gdB', itfType, noiseType, SNRdB), ...
                     'Position', [100 60 1200 980]);

        tiledlayout(4,1,'Padding','compact','TileSpacing','compact');

        ax1 = nexttile; hold(ax1,'on'); grid(ax1,'on');
        title(ax1, sprintf('DSSS: Variación H norm vs JSR | itf=%s | Ruido=%s | SNR=%gdB', upper(itfType), noiseType, SNRdB));
        xlabel(ax1,'JSR (dB)'); ylabel(ax1,'Delta H_norm');

        ax2 = nexttile; hold(ax2,'on'); grid(ax2,'on');
        title(ax2,'DSSS: KL(p_limpia || p_sucia) vs JSR');
        xlabel(ax2,'JSR (dB)'); ylabel(ax2,'KL (bits)');

        ax3 = nexttile; hold(ax3,'on'); grid(ax3,'on');
        title(ax3,'DSSS: MI I(X;Y) vs JSR');
        xlabel(ax3,'JSR (dB)'); ylabel(ax3,'MI (bits)');

        ax4 = nexttile; hold(ax4,'on'); grid(ax4,'on');
        title(ax4,'DSSS: Curtosis kappa_p vs JSR (sobre p_i normalizada)');
        xlabel(ax4,'JSR (dB)'); ylabel(ax4,'Curtosis');

        leg = strings(1,numel(f0_list));

        for ff = 1:numel(f0_list)
            f0 = f0_list(ff);
            cfgSig.f0 = f0;

            % ---- DSSS limpia (Ps=1) ----
            [~, x] = sig.generate('dsss', cfgSig);
            x = prep_unit_power(x);

            % ---- ruido elegido: baseline x+w ----
            cfgNoise = struct();
            cfgNoise.Fs     = Fs;
            cfgNoise.T      = T;
            cfgNoise.SNRdB  = SNRdB;
            cfgNoise.useDSP = true;
            cfgNoise.seed   = 1000 + 100*it + 10*nn + ff;

            try
                [~, w] = noi.generate(noiseType, cfgNoise);
            catch ME
                warning("No se pudo generar ruido '%s': %s", noiseType, ME.message);
                close(fig);
                continue;
            end

            w = w(:).';

            xr = x + w;         % referencia = señal+ruido
            xr = xr - mean(xr);

            % pmf espectral baseline
            p_base = psd_prob(xr, Fs, welch);
            H_base = shannon_H_norm(p_base);

            % para MI: estandarizar respecto al baseline
            mu0 = mean(xr);
            s0  = std(xr);
            if s0 == 0, s0 = 1; end

            % vectores de curvas
            dH = zeros(size(JSRdB_list));
            KL = zeros(size(JSRdB_list));
            MI = zeros(size(JSRdB_list));
            Kp = zeros(size(JSRdB_list));

            for jj = 1:numel(JSRdB_list)
                JSRdB = JSRdB_list(jj);

                cfgItf.JSRdB = JSRdB;
                cfgItf.seed  = 2000 + 100*it + 10*nn + 100*ff + jj;

                [~, jamm] = itf.generate(itfType, cfgItf);
                jamm = jamm(:).';

                y = xr + jamm;      % interferida
                y = y - mean(y);

                % ---- espectral ----
                p_y = psd_prob(y, Fs, welch);
                H_y = shannon_H_norm(p_y);

                dH(jj) = H_y - H_base;
                KL(jj) = kl_div(p_base, p_y);
                Kp(jj) = spectral_kurtosis(p_y);

                % ---- MI tiempo ----
                zr = (xr - mu0)/s0;
                zy = (y  - mu0)/s0;
                MI(jj) = mi_hist2(zr, zy, miCfg);
            end

            % ploteo
            plot(ax1, JSRdB_list, dH, 'LineWidth', 1.4);
            plot(ax2, JSRdB_list, KL, 'LineWidth', 1.4);
            plot(ax3, JSRdB_list, MI, 'LineWidth', 1.4);
            plot(ax4, JSRdB_list, Kp, 'LineWidth', 1.4);

            leg(ff) = sprintf('f0=%.0f Hz (%.3f Fs)', f0, f0/Fs);
        end

        legend(ax1, leg, 'Location','best');

        % Guardar
        fname = fullfile(outDir, sprintf('sec2_DSSS_%s_%s_SNR_%gdB.png', itfType, noiseType, SNRdB));
        exportgraphics(fig, fname);
        close(fig);

        fprintf("OK -> %s\n", fname);
    end
end

fprintf("\nListo. Figuras en: %s\n", outDir);

%% ===================== FUNCIONES LOCALES =====================

function x = prep_unit_power(x)
x = x(:).';
x = x - mean(x);
P = mean(x.^2);
if P == 0, P = 1; end
x = x / sqrt(P);
end

function p = psd_prob(x, Fs, welch)
x = x(:);
winLen  = min(welch.winLen, numel(x));
overlap = min(welch.overlap, winLen-1);
nfft    = max(welch.nfft, 2^nextpow2(winLen));
[Pxx, ~] = pwelch(x, winLen, overlap, nfft, Fs, 'onesided');
Pxx = max(Pxx, 0);
p = Pxx / sum(Pxx);
epsVal = 1e-12;
p = max(p, epsVal);
p = p / sum(p);
end

function Hn = shannon_H_norm(p)
p = p(:);
H  = -sum(p .* log2(p));
M  = numel(p);
Hn = H / log2(M);
end

function D = kl_div(p, q)
p = p(:); q = q(:);
epsVal = 1e-12;
p = max(p, epsVal); p = p/sum(p);
q = max(q, epsVal); q = q/sum(q);
D = sum(p .* log2(p ./ q));
end

function k = spectral_kurtosis(p)
p = p(:);
M  = numel(p);
mu = 1/M;
v = mean((p - mu).^2);
if v <= 1e-30
    k = NaN;
    return;
end
k = mean((p - mu).^4) / (v^2);
end

function I = mi_hist2(x, y, cfg)
x = x(:); y = y(:);
L = cfg.L;
B = cfg.B;
x = min(max(x, -L), L);
y = min(max(y, -L), L);
edges = linspace(-L, L, B+1);
Nxy = histcounts2(x, y, edges, edges);
Nxy = Nxy + cfg.alpha;
pxy = Nxy / sum(Nxy(:));
px  = sum(pxy, 2);
py  = sum(pxy, 1);
den  = px * py;
mask = pxy > 0 & den > 0;
I = sum(pxy(mask) .* log2(pxy(mask) ./ den(mask)));
end


