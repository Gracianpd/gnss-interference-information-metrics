function out = ti_interference()
% TI_INTERFERENCE  Generación y plot de interferencias (CW, NBI, chirp).
%
% Uso:
%   itf = ti_interference();
%   cfg = struct('Fs',10000,'T',1.0,'JSRdB',0,'seed',1,'useDSP',true, ...
%                'fj',2000,'Bj',400,'f1',500,'f2',4500,'phi',0);
%   J = itf.generate_all(cfg);
%   itf.plot_all(J);

out.generate        = @generate;
out.generate_all    = @generate_all;
out.plot_all        = @plot_all;
out.scale_to_power  = @scale_to_power;
end

% ----------------- Utils -----------------

function j = scale_to_power(j, Pj)
% Quita media y escala a potencia media discreta Pj
j = j(:);
j = j - mean(j);
Pj0 = mean(j.^2);
j = j * sqrt(Pj / Pj0);
end

% ----------------- Generación individual -----------------

function [t, j] = generate(type, cfg)
% GENERATE  Genera interferencia discreta j[n].

% Defaults
if ~isfield(cfg,'useDSP'), cfg.useDSP = true; end
if ~isfield(cfg,'seed'), cfg.seed = 1; end
if ~isfield(cfg,'phi'),  cfg.phi = 0; end

Fs = cfg.Fs;
T  = cfg.T;
N  = round(Fs*T);
t  = (0:N-1)/Fs;

Pj = 10^(cfg.JSRdB/10);  % asumiendo Ps=1

type = lower(strtrim(type));

% DSP availability
hasBandpass = (exist('dsp.BandpassFilter','class') == 8);
hasChirpObj = (exist('dsp.Chirp','class') == 8);

rng(cfg.seed);

switch type
    case {'cw'}
        if ~isfield(cfg,'fj'), error("CW requiere cfg.fj (Hz)."); end
        fj = cfg.fj;
        j0 = cos(2*pi*fj*t + cfg.phi);
        j  = scale_to_power(j0, Pj);

    case {'nbi','narrowband'}
        if ~isfield(cfg,'fj') || ~isfield(cfg,'Bj')
            error("NBI requiere cfg.fj (Hz) y cfg.Bj (Hz).");
        end
        fj = cfg.fj;
        Bj = cfg.Bj;

        % Validación simple de banda
        if fj - Bj/2 <= 0 || fj + Bj/2 >= Fs/2
            error("Elige fj y Bj para que (fj-Bj/2)>0 y (fj+Bj/2)<Fs/2.");
        end

        u = randn(N,1);

        if cfg.useDSP && hasBandpass
            % Pasabanda DSP (parámetros conservadores)
            bp = dsp.BandpassFilter( ...
                'SampleRate', Fs, ...
                'StopbandFrequency1', max(1, fj - Bj), ...
                'PassbandFrequency1', fj - Bj/2, ...
                'PassbandFrequency2', fj + Bj/2, ...
                'StopbandFrequency2', min(Fs/2-1, fj + Bj));
            j0 = bp(u);
        else
            % Fallback: designfilt + filter
            d = designfilt('bandpassiir', ...
                'FilterOrder', 6, ...
                'HalfPowerFrequency1', fj - Bj/2, ...
                'HalfPowerFrequency2', fj + Bj/2, ...
                'SampleRate', Fs);
            j0 = filter(d, u);
        end

        j  = scale_to_power(j0, Pj);

    case {'chirp'}
        if ~isfield(cfg,'f1') || ~isfield(cfg,'f2')
            error("Chirp requiere cfg.f1 y cfg.f2 (Hz).");
        end
        f1 = cfg.f1;
        f2 = cfg.f2;

        if f1 < 0 || f2 < 0 || f1 >= Fs/2 || f2 >= Fs/2
            error("f1 y f2 deben estar en [0, Fs/2).");
        end

        if cfg.useDSP && hasChirpObj
            ch = dsp.Chirp( ...
                Type="Linear", ...
                SweepDirection="Unidirectional", ...
                InitialFrequency=f1, ...
                TargetFrequency=f2, ...
                TargetTime=T, ...
                SweepTime=T, ...
                InitialPhase=cfg.phi, ...
                SamplesPerFrame=N, ...
                SampleRate=Fs);
            j0 = ch();   % Nx1
        else
            % Fallback por fórmula: fase de chirp lineal
            % f(t) = f1 + k t, k = (f2-f1)/T
            tt = t(:);
            k  = (f2 - f1)/T;
            theta = 2*pi*(f1*tt + 0.5*k*tt.^2) + cfg.phi;
            j0 = cos(theta);
        end

        j = scale_to_power(j0, Pj);

    otherwise
        error("Tipo de interferencia '%s' no reconocido.", type);
end

j = j.'; % fila
end

% ----------------- Generación de todas -----------------

function J = generate_all(cfg)
% GENERATE_ALL  Genera CW, NBI y chirp

types = {'cw','nbi','chirp'};
names = {'cw','nbi','chirp'};

J = struct();

[t, jcw] = generate('cw', cfg);
J.t  = t;
J.Fs = cfg.Fs;
J.JSRdB = cfg.JSRdB;

J.(names{1}) = jcw;

for k = 2:numel(types)
    [~, jk] = generate(types{k}, cfg);
    J.(names{k}) = jk;
end

% guardar parámetros útiles para títulos
J.fj = cfg.fj;
if isfield(cfg,'Bj'), J.Bj = cfg.Bj; end
if isfield(cfg,'f1'), J.f1 = cfg.f1; end
if isfield(cfg,'f2'), J.f2 = cfg.f2; end

end

% ----------------- Plot conjunto -----------------

function fig = plot_all(J)
% PLOT_ALL  Plotea interferencias en tiempo y PSD (Welch) en lineal.

t  = J.t;
Fs = J.Fs;

fields = {'cw','nbi','chirp'};
titles = { ...
    sprintf('CW (f_j=%.1f Hz)', J.fj), ...
    sprintf('NBI (f_j=%.1f Hz, B=%.1f Hz)', J.fj, J.Bj), ...
    sprintf('Chirp (f_1=%.1f Hz -> f_2=%.1f Hz)', J.f1, J.f2) ...
};

% Para tiempo: mostrar unas muestras (p.ej. primeros 3000)
nShow = min(numel(t), 3000);

% Welch (lineal)
N = numel(t);
winLen  = min(4096, 2^floor(log2(max(256, floor(N/4)))));
overlap = floor(0.5*winLen);
nfft    = max(4096, 2^nextpow2(winLen));

fig = figure('Name','Interferencias','Position',[100 100 1200 800]);

% Tiempo
for i = 1:3
    subplot(3,2,2*i-1);
    j = J.(fields{i});
    plot(t(1:nShow), j(1:nShow), 'LineWidth', 1);
    grid on;
    title(sprintf('%s | JSRdB=%.1f', titles{i}, J.JSRdB));
    xlabel('t (s)'); ylabel('j[n]');
end

% PSD (Welch) - LINEAL (sin dB)
for i = 1:3
    subplot(3,2,2*i);
    j = J.(fields{i});
    [Pxx, f] = pwelch(j(:), winLen, overlap, nfft, Fs, 'onesided');
    plot(f, Pxx, 'LineWidth', 1);
    grid on;
    title(['PSD Welch - ' titles{i}]);
    xlabel('f (Hz)'); ylabel('PSD');
end

end