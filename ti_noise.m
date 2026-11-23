function out = ti_noise()
% TI_NOISE  Generación y plot de ruidos (white/pink/brown).
%
% Uso:
%   n = ti_noise();
%   cfg = struct('Fs',10000,'T',1.0,'SNRdB',10,'seed',1,'useDSP',true);
%   W = n.generate_all(cfg);
%   n.plot_all(W);

out.generate      = @generate;
out.generate_all  = @generate_all;
out.plot_all      = @plot_all;
out.scale_to_power = @scale_to_power;
end

% ----------------- Utils -----------------

function w = scale_to_power(w, Pw)
% Quita media y escala a potencia media discreta Pw
w = w(:);
w = w - mean(w);
w = w * sqrt(Pw / mean(w.^2));
end

% ----------------- Generación individual -----------------

function [t, w] = generate(type, cfg)
% GENERATE  Genera ruido discreto w[n] del tipo indicado.

% Defaults
if ~isfield(cfg,'useDSP'), cfg.useDSP = true; end
if ~isfield(cfg,'seed'), cfg.seed = 1; end

Fs  = cfg.Fs;
T   = cfg.T;
N   = round(Fs*T);
t   = (0:N-1)/Fs;

SNRdB = cfg.SNRdB;
Pw = 10^(-SNRdB/10);  % asumiendo señal normalizada a potencia 1

rng(cfg.seed);

type = lower(strtrim(type));
hasColoredNoise = (exist('dsp.ColoredNoise','class') == 8);

switch type
    case {'white','awgn'}
        w = randn(N,1);
        w = scale_to_power(w, Pw);

    case {'pink'}
        if cfg.useDSP && hasColoredNoise
            cn = dsp.ColoredNoise('Color','pink','SamplesPerFrame',N,'NumChannels',1);
            w = cn();
            w = scale_to_power(w, Pw);
        else
            error("Pink noise requiere DSP System Toolbox (dsp.ColoredNoise) o implementar fallback.");
        end

    case {'brown','brownian'}
        if cfg.useDSP && hasColoredNoise
            cn = dsp.ColoredNoise('Color','brown','SamplesPerFrame',N,'NumChannels',1);
            w = cn();
            w = scale_to_power(w, Pw);
        else
            error("Brown noise requiere DSP System Toolbox (dsp.ColoredNoise) o implementar fallback.");
        end

    otherwise
        error("Tipo de ruido '%s' no reconocido.", type);
end

w = w.'; % devolver fila para consistencia con tu estilo
end

% ----------------- Generación de todos -----------------

function W = generate_all(cfg)
% GENERATE_ALL  Genera white/pink/brown y devuelve struct con t y ruidos.

noiseTypes = {'white','pink','brown'};
names      = {'white','pink','brown'};

W = struct();

[t, w0] = generate('white', cfg);
W.t  = t;
W.Fs = cfg.Fs;
W.SNRdB = cfg.SNRdB;
W.white = w0;

for k = 2:numel(noiseTypes)
    [~, wk] = generate(noiseTypes{k}, cfg);
    W.(names{k}) = wk;
end

end

% ----------------- Plot conjunto -----------------

function fig = plot_all(W)
% PLOT_ALL  Plotea en tiempo y PSD (Welch) para white/pink/brown.

t = W.t;
Fs = W.Fs;

fields = {'white','pink','brown'};
titles = {'Ruido blanco','Ruido rosa (~1/f)','Ruido marrón (~1/f^2)'};

% Mostrar solo una parte en tiempo
nShow = min(numel(t), 3000);

fig = figure('Name','Ruidos','Position',[100 100 1200 800]);

% ---- Tiempo
for i = 1:3
    subplot(3,2,2*i-1);
    w = W.(fields{i});
    plot(t(1:nShow), w(1:nShow), 'LineWidth', 1);
    grid on;
    title(sprintf('%s | SNRdB=%.1f', titles{i}, W.SNRdB));
    xlabel('t (s)'); ylabel('w[n]');
end

% ---- PSD (Welch)
N = numel(t);
winLen  = min(4096, 2^floor(log2(max(256, floor(N/4)))));
overlap = floor(0.5*winLen);
nfft    = max(4096, 2^nextpow2(winLen));

for i = 1:3
    subplot(3,2,2*i);
    w = W.(fields{i});
    [Pxx, f] = pwelch(w(:), winLen, overlap, nfft, Fs, 'onesided');
    plot(f, 10*log10(Pxx), 'LineWidth', 1);
    grid on;
    title(['PSD Welch - ' titles{i}]);
    xlabel('f (Hz)'); ylabel('PSD (dB/Hz)');
end

end