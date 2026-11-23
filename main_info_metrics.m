%% main_info_metrics.m
% Calcula Shannon (espectral), KL (espectral) y MI clásica (tiempo)
% para todas las combinaciones: señal x ruido x interferencia.

clear; close all; clc;

%% ====== Instanciar módulos ======
sig = ti_signals();

noi = ti_noise();

itf = ti_interference();

%% ====== Config global ======
cfgBase = struct();
cfgBase.Fs  = 10000;
cfgBase.T   = 1.0;          % 1 s para PSD y MI más estables
cfgBase.A   = 1;
cfgBase.phi = 0;

% Señales
cfgBase.f0   = 0.02*cfgBase.Fs;   % frecuencia base de señal
cfgBase.D    = 0.25;              % duty rectangular
cfgBase.Rc   = cfgBase.Fs/10;     % DSSS chip-rate (Nc=10)
cfgBase.seed = 1;

% Ruido
cfgNoise = struct();
cfgNoise.Fs = cfgBase.Fs;
cfgNoise.T  = cfgBase.T;
cfgNoise.SNRdB = 10;          % <-- ajusta aquí (o haz lista)
cfgNoise.useDSP = true;

% Interferencias
cfgItf = struct();
cfgItf.Fs = cfgBase.Fs;
cfgItf.T  = cfgBase.T;
cfgItf.JSRdB = 0;             % <-- ajusta aquí (o haz lista)
cfgItf.useDSP = true;
cfgItf.phi = 0;

% CW/NBI
cfgItf.fj = 2000;
cfgItf.Bj = 400;

% Chirp
cfgItf.f1 = 500;
cfgItf.f2 = 4500;

% Welch (para tu p_i)
welch = struct();
welch.winLen  = 4096;
welch.overlap = 2048;
welch.nfft    = 4096;

% MI (histograma 2D)
miCfg = struct();
miCfg.B = 32;         % bins por eje
miCfg.alpha = 1e-6;   % pseudoconteo (suavizado)
miCfg.pLo = 1;        % percentil low para edges
miCfg.pHi = 99;       % percentil high para edges

%% ====== Listas de combinaciones ======
signalTypes = {'sin','sq','rect','tri','saw','dsss'};
noiseTypes  = {'none','white','pink','brown'};
itfTypes    = {'none','cw','nbi','chirp'};

%% ====== Bucle ======
rows = [];
rowIdx = 0;

for s = 1:numel(signalTypes)
    sigType = signalTypes{s};

    % Generar señal limpia
    [t, x] = sig.generate(sigType, cfgBase);
    x = prep_signal(x);              % quitar media + potencia unitaria (Ps=1)

    % Distribución espectral de referencia (limpia)
    p_clean = psd_prob(x, cfgBase.Fs, welch);
    H_clean = shannon_H(p_clean);

    for n = 1:numel(noiseTypes)
        nType = noiseTypes{n};

        % Ruido
        w = zeros(size(x));
        if ~strcmp(nType,'none')
            cfgNoise.seed = 1000 + s*100 + n;  % reproducible por combinación
            try
                [~, w] = noi.generate(nType, cfgNoise);
            catch ME
                warning("Saltando ruido '%s' (prob. toolbox/licencia): %s", nType, ME.message);
                continue;
            end
            w = w(:).'; % asegurar fila
        end

        for jidx = 1:numel(itfTypes)
            jType = itfTypes{jidx};

            % Interferencia
            j = zeros(size(x));
            if ~strcmp(jType,'none')
                cfgItf.seed = 2000 + s*100 + jidx;
                try
                    [~, j] = itf.generate(jType, cfgItf);
                catch ME
                    warning("Saltando interferencia '%s': %s", jType, ME.message);
                    continue;
                end
                j = j(:).';
            end

            % Señal recibida
            y = x + w + j;
            %y = prep_signal(y);   % coherente con tu pipeline (quita DC, normaliza)
            y = y- mean(y)
            % ===== Métricas =====
            % Shannon (espectral) de y
            p_y = psd_prob(y, cfgBase.Fs, welch);
            H_y = shannon_H(p_y);

            % KL (espectral) limpia -> afectada
            Dkl = kl_div(p_clean, p_y);

            % MI clásica (tiempo) entre x e y
            Ixy = mi_hist2(x, y, miCfg);

            % Guardar fila
            rowIdx = rowIdx + 1;
            rows(rowIdx).signal = sigType;
            rows(rowIdx).noise  = nType;
            rows(rowIdx).itf    = jType;
            rows(rowIdx).Fs     = cfgBase.Fs;
            rows(rowIdx).T      = cfgBase.T;
            rows(rowIdx).f0     = cfgBase.f0;
            rows(rowIdx).SNRdB  = cfgNoise.SNRdB;
            rows(rowIdx).JSRdB  = cfgItf.JSRdB;
            rows(rowIdx).H_clean = H_clean;
            rows(rowIdx).H_y     = H_y;
            rows(rowIdx).KL_pq   = Dkl;
            rows(rowIdx).MI_xy   = Ixy;
        end
    end
end

% Pasar a tabla y mostrar
T = struct2table(rows);
disp(T);

% Guardar resultados
if ~exist('Results','dir'), mkdir('Results'); end
writetable(T, fullfile('Results','info_metrics.csv'));
save(fullfile('Results','info_metrics.mat'), 'T', 'cfgBase', 'cfgNoise', 'cfgItf', 'welch', 'miCfg');

fprintf("\nListo. Guardado en Results/info_metrics.csv y Results/info_metrics.mat\n");

%% ====== Funciones locales ======

function z = prep_signal(z)
% Quitar media y normalizar potencia a 1 (Ps=1)
z = z(:).';
z = z - mean(z);
P = mean(z.^2);
if P == 0, P = 1; end
z = z / sqrt(P);
end

function p = psd_prob(x, Fs, welch)
% PSD Welch -> pmf espectral p_i = P_i / sum(P_i)
x = x(:);
[Pxx, ~] = pwelch(x, welch.winLen, welch.overlap, welch.nfft, Fs, 'onesided');
Pxx = max(Pxx, 0);
p = Pxx / sum(Pxx);
% regularización por seguridad numérica
epsVal = 1e-12;
p = max(p, epsVal);
p = p / sum(p);
end

function H = shannon_H(p)
% Entropía de Shannon en bits
p = p(:);
H = -sum(p .* log2(p));
end

function D = kl_div(p, q)
% KL(p||q) en bits con regularización
p = p(:); q = q(:);
epsVal = 1e-12;
p = max(p, epsVal); p = p/sum(p);
q = max(q, epsVal); q = q/sum(q);
D = sum(p .* log2(p ./ q));
end

function I = mi_hist2(x, y, miCfg)
% MI clásica I(X;Y) estimada por histograma 2D (bits)
x = x(:); y = y(:);

% estandarizar (robustez numérica)
%x = (x - mean(x)) / std(x);
%y = (y - mean(y)) / std(y);
x = x - mean(x);
y = y - mean(y);


% edges por percentiles para evitar outliers
loX = prctile(x, miCfg.pLo); hiX = prctile(x, miCfg.pHi);
loY = prctile(y, miCfg.pLo); hiY = prctile(y, miCfg.pHi);

% evitar degeneración si std~0
if hiX == loX, hiX = loX + 1; end
if hiY == loY, hiY = loY + 1; end

edgesX = linspace(loX, hiX, miCfg.B+1);
edgesY = linspace(loY, hiY, miCfg.B+1);

% recortar para que todo caiga en bins
x = min(max(x, loX), hiX);
y = min(max(y, loY), hiY);

Nxy = histcounts2(x, y, edgesX, edgesY);

% suavizado (pseudoconteo)
Nxy = Nxy + miCfg.alpha;

pxy = Nxy / sum(Nxy(:));
px  = sum(pxy, 2);
py  = sum(pxy, 1);

den = px * py;               % BxB
mask = pxy > 0 & den > 0;

I = sum( pxy(mask) .* log2( pxy(mask) ./ den(mask) ) );
end


