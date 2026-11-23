%% main_sweep_metrics_and_plots.m
% Barrido SNR/JSR + métricas + (1) curvas vs JSR (2) heatmaps SNR-JSR (3) ROC/AUC

clear; close all; clc;

%% ====== Módulos ======
sig = ti_signals();
noi = ti_noise();
itf = ti_interference();

%% ====== Config base ======
cfgBase = struct();
cfgBase.Fs  = 10000;
cfgBase.T   = 1.0;
cfgBase.A   = 1;
cfgBase.phi = 0;
cfgBase.f0   = 0.02*cfgBase.Fs;
cfgBase.D    = 0.25;
cfgBase.Rc   = cfgBase.Fs/10;
cfgBase.seed = 1;

cfgNoise = struct('Fs',cfgBase.Fs,'T',cfgBase.T,'useDSP',true,'seed',1,'SNRdB',10);
cfgItf   = struct('Fs',cfgBase.Fs,'T',cfgBase.T,'useDSP',true,'seed',1,'JSRdB',0,'phi',0, ...
                  'fj',2000,'Bj',400,'f1',500,'f2',4500);

welch = struct('winLen',4096,'overlap',2048,'nfft',4096);

miCfg = struct();
miCfg.B = 32;
miCfg.alpha = 1e-6;
miCfg.pLo = 1;
miCfg.pHi = 99;

%% ====== Barridos ======
SNRdB_list = [0 5 10 15 20];          % ajusta
JSRdB_list = [-30 -20 -10 0 10 20];   % ajusta

signalTypes = {'sin','sq','rect','tri','saw','dsss'};
noiseTypes  = {'none','white','pink','brown'};
itfTypes    = {'none','cw','nbi','chirp'};

%% ====== Run ======
rows = [];
rowIdx = 0;

for s = 1:numel(signalTypes)
    sigType = signalTypes{s};

    % Señal limpia x (normalizada a potencia 1)
    [t, xraw] = sig.generate(sigType, cfgBase);
    x = normalize_unit_power(xraw);

    % Referencias espectrales de x
    p_clean = psd_prob(x, cfgBase.Fs, welch);
    Hc = shannon_H_norm(p_clean);

    for n = 1:numel(noiseTypes)
        nType = noiseTypes{n};

        if strcmp(nType,'none')
            SNRs = NaN;   % no barrer, no aplica
        else
            SNRs = SNRdB_list;
        end

        for sn = 1:numel(SNRs)
            cfgNoise.SNRdB = SNRs(sn);
            cfgNoise.seed = 1000 + s*100 + n*10 + sn;

            % Ruido
            w = zeros(size(x));
            if ~strcmp(nType,'none')
                [~, w] = noi.generate(nType, cfgNoise);
                w = w(:).';
            end

            for jidx = 1:numel(itfTypes)
                jType = itfTypes{jidx};

                if strcmp(jType,'none')
                    JSRs = NaN;  % no barrer, no aplica
                else
                    JSRs = JSRdB_list;
                end

                for jj = 1:numel(JSRs)
                    cfgItf.JSRdB = JSRs(jj);
                    cfgItf.seed = 2000 + s*100 + jidx*10 + jj;

                    % Interferencia
                    j = zeros(size(x));
                    if ~strcmp(jType,'none')
                        [~, j] = itf.generate(jType, cfgItf);
                        j = j(:).';
                    end

                    % Señal recibida (NO normalizar potencia: mantiene SNR/JSR)
                    y = x + w + j;
                    y = y - mean(y);

                    % ===== Métricas =====
                    p_y  = psd_prob(y, cfgBase.Fs, welch);
                    Hy   = shannon_H_norm(p_y);
                    dH   = Hy - Hc;

                    Dkl  = kl_div(p_clean, p_y);

                    Ixy  = mi_hist2_fixed_edges(x, y, miCfg);

                    % Guardar
                    rowIdx = rowIdx + 1;
                    rows(rowIdx).signal  = string(sigType);
                    rows(rowIdx).noise   = string(nType);
                    rows(rowIdx).itf     = string(jType);
                    rows(rowIdx).Fs      = cfgBase.Fs;
                    rows(rowIdx).T       = cfgBase.T;
                    rows(rowIdx).f0      = cfgBase.f0;
                    rows(rowIdx).SNRdB   = cfgNoise.SNRdB;
                    rows(rowIdx).JSRdB   = cfgItf.JSRdB;

                    rows(rowIdx).H_clean = Hc;
                    rows(rowIdx).H_y     = Hy;
                    rows(rowIdx).dH      = dH;

                    rows(rowIdx).KL_pq   = Dkl;
                    rows(rowIdx).MI_xy   = Ixy;
                end
            end
        end
    end
end

T = struct2table(rows);
disp(T(1:min(height(T),10),:));

if ~exist('Results','dir'), mkdir('Results'); end
writetable(T, fullfile('Results','info_metrics_sweep.csv'));
save(fullfile('Results','info_metrics_sweep.mat'), 'T', 'cfgBase', 'cfgNoise', 'cfgItf', 'welch', 'miCfg');

%% ====== PLOTS (las 3 cosas) ======
if ~exist('Figures','dir'), mkdir('Figures'); end

% (1) Curvas vs JSR (sensibilidad): ejemplo con 3 señales y ruido fijo
plotCurvesCfg = struct();
plotCurvesCfg.signals = ["sin","sq","dsss"];
plotCurvesCfg.noise   = "white";     % cambia a "none"/"pink"/...
plotCurvesCfg.SNRdB   = 10;          % elige un SNR concreto
plotCurvesCfg.itfs    = ["cw","nbi","chirp"];
plot_curves_vs_jsr(T, plotCurvesCfg);

% (2) Heatmaps SNR–JSR para DSSS
plotHeatCfg = struct();
plotHeatCfg.signal = "dsss";
plotHeatCfg.noises = ["white","pink","brown"];   % puedes dejar solo "white"
plotHeatCfg.itfs   = ["cw","nbi","chirp"];
plot_heatmaps_snr_jsr(T, plotHeatCfg);

% (3) ROC + AUC global (itf vs no-itf)
plot_roc_auc(T);

fprintf("\nOK: Results/info_metrics_sweep.csv + Figures/*.png\n");

%% ===================== FUNCIONES LOCALES =====================

function x = normalize_unit_power(x)
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
H = -sum(p .* log2(p));
M = numel(p);
Hn = H / log2(M);
end

function D = kl_div(p, q)
p = p(:); q = q(:);
epsVal = 1e-12;
p = max(p, epsVal); p = p/sum(p);
q = max(q, epsVal); q = q/sum(q);
D = sum(p .* log2(p ./ q));
end

function I = mi_hist2_fixed_edges(x, y, miCfg)
% MI con edges FIJOS definidos por x (referencia), para no "adaptar" bins a y
x = x(:); y = y(:);
x = x - mean(x);
y = y - mean(y);

lo = prctile(x, miCfg.pLo);
hi = prctile(x, miCfg.pHi);
if hi == lo, hi = lo + 1; end

edges = linspace(lo, hi, miCfg.B+1);

% recortar a rango
x = min(max(x, lo), hi);
y = min(max(y, lo), hi);

Nxy = histcounts2(x, y, edges, edges);
Nxy = Nxy + miCfg.alpha;

pxy = Nxy / sum(Nxy(:));
px  = sum(pxy, 2);
py  = sum(pxy, 1);

den = px * py;
mask = pxy > 0 & den > 0;
I = sum( pxy(mask) .* log2( pxy(mask) ./ den(mask) ) );
end

%% ---------------- (1) Curvas vs JSR ----------------
function plot_curves_vs_jsr(T, cfg)
% Curvas de KL, MI y dH vs JSR para varias interferencias, señal y ruido fijos
sigList = cfg.signals;
noise   = cfg.noise;
SNRdB   = cfg.SNRdB;
itfs    = cfg.itfs;

for s = 1:numel(sigList)
    sig = sigList(s);

    % Filtrar por señal + ruido + SNR
    TT = T(T.signal==sig & T.noise==noise & isfinite(T.SNRdB) & abs(T.SNRdB - SNRdB)<1e-9, :);
    if isempty(TT), warning("No hay datos para %s / %s / SNR=%.1f", sig, noise, SNRdB); continue; end

    % Figura 3 subplots
    fig = figure('Name',sprintf('Curvas vs JSR - %s', sig),'Position',[100 100 1100 700]);
    tiledlayout(3,1,'Padding','compact');

    % KL
    nexttile; hold on; grid on;
    for k = 1:numel(itfs)
        itf = itfs(k);
        S = TT(TT.itf==itf & isfinite(TT.JSRdB), :);
        [jsr, idx] = sort(S.JSRdB);
        plot(jsr, S.KL_pq(idx), 'LineWidth', 1.5, 'DisplayName', char(itf));
    end
    xlabel('JSR (dB)'); ylabel('KL(p||q) [bits]');
    title(sprintf('%s | noise=%s | SNR=%.1f dB | KL', sig, noise, SNRdB));
    legend('Location','best');

    % MI
    nexttile; hold on; grid on;
    for k = 1:numel(itfs)
        itf = itfs(k);
        S = TT(TT.itf==itf & isfinite(TT.JSRdB), :);
        [jsr, idx] = sort(S.JSRdB);
        plot(jsr, S.MI_xy(idx), 'LineWidth', 1.5, 'DisplayName', char(itf));
    end
    xlabel('JSR (dB)'); ylabel('I(X;Y) [bits]');
    title('Mutual Information (tiempo)');
    legend('Location','best');

    % dH
    nexttile; hold on; grid on;
    for k = 1:numel(itfs)
        itf = itfs(k);
        S = TT(TT.itf==itf & isfinite(TT.JSRdB), :);
        [jsr, idx] = sort(S.JSRdB);
        plot(jsr, S.dH(idx), 'LineWidth', 1.5, 'DisplayName', char(itf));
    end
    xlabel('JSR (dB)'); ylabel('\Delta \hat{H} (norm.)');
    title('Cambio de entropía espectral normalizada');
    legend('Location','best');

    exportgraphics(fig, fullfile('Figures',sprintf('curves_vs_jsr_%s_noise_%s_SNR_%g.png',sig,noise,SNRdB)));
end
end

%% ---------------- (2) Heatmaps SNR–JSR ----------------
function plot_heatmaps_snr_jsr(T, cfg)
sig   = cfg.signal;
noises= cfg.noises;
itfs  = cfg.itfs;

SNRs = unique(T.SNRdB(isfinite(T.SNRdB)));
JSRs = unique(T.JSRdB(isfinite(T.JSRdB)));

for nn = 1:numel(noises)
    noise = noises(nn);
    for ii = 1:numel(itfs)
        itf = itfs(ii);

        TT = T(T.signal==sig & T.noise==noise & T.itf==itf & isfinite(T.SNRdB) & isfinite(T.JSRdB), :);
        if isempty(TT), warning("No hay datos heatmap para %s/%s/%s",sig,noise,itf); continue; end

        % Construir matrices
        Z_KL = nan(numel(JSRs), numel(SNRs));
        Z_MI = nan(numel(JSRs), numel(SNRs));
        Z_dH = nan(numel(JSRs), numel(SNRs));

        for a = 1:numel(JSRs)
            for b = 1:numel(SNRs)
                idx = TT.JSRdB==JSRs(a) & TT.SNRdB==SNRs(b);
                if any(idx)
                    % si hay duplicados, usar media
                    Z_KL(a,b) = mean(TT.KL_pq(idx));
                    Z_MI(a,b) = mean(TT.MI_xy(idx));
                    Z_dH(a,b) = mean(TT.dH(idx));
                end
            end
        end

        fig = figure('Name',sprintf('Heatmaps %s %s %s',sig,noise,itf),'Position',[100 100 1200 800]);
        tiledlayout(1,3,'Padding','compact');

        nexttile; imagesc(SNRs, JSRs, Z_KL); axis xy; colorbar;
        xlabel('SNR (dB)'); ylabel('JSR (dB)'); title('KL(p||q)');

        nexttile; imagesc(SNRs, JSRs, Z_MI); axis xy; colorbar;
        xlabel('SNR (dB)'); ylabel('JSR (dB)'); title('MI I(X;Y)');

        nexttile; imagesc(SNRs, JSRs, Z_dH); axis xy; colorbar;
        xlabel('SNR (dB)'); ylabel('JSR (dB)'); title('\Delta \hat{H}');

        exportgraphics(fig, fullfile('Figures',sprintf('heatmap_%s_noise_%s_itf_%s.png',sig,noise,itf)));
    end
end
end

%% ---------------- (3) ROC + AUC ----------------
function plot_roc_auc(T)
% Clasificación: interferencia presente vs no
isPos = T.itf ~= "none";     % positivos: hay interferencia
isNeg = T.itf == "none";     % negativos: no hay interferencia

% Para ROC, elimina filas NaN (por ejemplo combinaciones sin barrido)
ok = ~isnan(T.KL_pq) & ~isnan(T.MI_xy) & ~isnan(T.dH);
T2 = T(ok,:);

labels = (T2.itf ~= "none"); % 1 = interferencia

% Scores (más alto => más probable interferencia)
scoreKL = T2.KL_pq;
scoreMI = -T2.MI_xy;         % MI baja => interferencia
scoredH = abs(T2.dH);        % cambio de entropía (módulo)

% ROC + AUC (con perfcurve si está disponible)
hasPerf = exist('perfcurve','file') == 2;

fig = figure('Name','ROC/AUC','Position',[100 100 900 700]);
hold on; grid on;

if hasPerf
    [X1,Y1,~,A1] = perfcurve(labels, scoreKL, true);
    [X2,Y2,~,A2] = perfcurve(labels, scoreMI, true);
    [X3,Y3,~,A3] = perfcurve(labels, scoredH, true);
else
    [X1,Y1,A1] = roc_manual(labels, scoreKL);
    [X2,Y2,A2] = roc_manual(labels, scoreMI);
    [X3,Y3,A3] = roc_manual(labels, scoredH);
end

plot(X1,Y1,'LineWidth',2,'DisplayName',sprintf('KL (AUC=%.3f)',A1));
plot(X2,Y2,'LineWidth',2,'DisplayName',sprintf('-MI (AUC=%.3f)',A2));
plot(X3,Y3,'LineWidth',2,'DisplayName',sprintf('|dH| (AUC=%.3f)',A3));

plot([0 1],[0 1],'k--','DisplayName','Azar');

xlabel('False Positive Rate'); ylabel('True Positive Rate');
title('ROC para detección de interferencia (itf \neq none)');
legend('Location','southeast');

exportgraphics(fig, fullfile('Figures','roc_auc.png'));
end

function [FPR,TPR,AUC] = roc_manual(labels, scores)
% ROC manual + AUC por trapecios (si no hay perfcurve)
labels = logical(labels(:));
scores = scores(:);

% ordenar por score desc
[sc, idx] = sort(scores,'descend');
lab = labels(idx);

P = sum(lab);
N = sum(~lab);

TP = cumsum(lab);
FP = cumsum(~lab);

TPR = [0; TP/P; 1];
FPR = [0; FP/N; 1];

% AUC trapecios
AUC = trapz(FPR, TPR);
end


