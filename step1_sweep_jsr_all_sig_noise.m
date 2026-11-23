function step1_sweep_jsr_all_sig_noise()
% STEP1_SWEEP_JSR_ALL_SIG_NOISE
% Barrido en JSR para 3 interferencias (cw,nbi,chirp), para todas las
% combinaciones (signal, noise) SIN caso "sin ruido".
%
% PLOTS (1 figura por signal+noise):
%   - dH = Hhat(PSD(x+w+j)) - Hhat(PSD(x+w))   (NO entropia absoluta)
%   - MI
%   - H(Y|X) (ruido del canal)
%   - Eff = I/H(Y)
%   - KL(p_x || p_y)
%
% CSV:
%   - csv_full_jsr_sweep.csv: todas filas + baseline solo ruido (columnas *_noiseonly)
%   - csv_rep_jsr.csv: resumen para JSR_rep + baseline solo ruido
%
% Todo se vuelca en carpeta: outputs_seccion1/
%
% Requisitos:
%   - ti_signals.m
%   - ti_noise.m
%   - ti_interference.m

clc; close all;

%% ===== Carpeta outputs =====
outDir = 'outputs_seccion1';
if ~exist(outDir,'dir')
    mkdir(outDir);
end

%% ===== Config global =====
cfgBase = struct();
cfgBase.Fs   = 10000;
cfgBase.T    = 1.0;
cfgBase.A    = 1;
cfgBase.phi  = 0;
cfgBase.D    = 0.25;     % rectangular
cfgBase.Rc   = 1000;     % DSSS (Fs/Rc = 10 entero)
cfgBase.seed = 1;

% Frecuencia fija (representativa)
f0 = 1000;
cfgBase.f0 = f0;

% Ruido fijo (para este bloque)
SNRdB = 3;

% Barrido JSR
JSRdB_list = -30:5:10;

% Interferencias
itfList = {'cw','nbi','chirp'};

% Parametros NBI / Chirp
Bj = 400;          % Hz
chirpDelta = 500;  % chirp alrededor de f0

% Welch / MI
Wcfg.winLen  = 2048;
Wcfg.overlap = 1024;
Wcfg.nfft    = 4096;
epsPSD = 1e-12;

MIcfg.L = 6;
MIcfg.B = 64;
MIcfg.eps = 1e-12;
MIcfg.Neff_max = 20000;

% Listas (sin 'none')
signals = {'sin','sq','rect','tri','saw','dsss'};
noises  = {'white','pink','brown'};

% JSR representativa para el CSV resumen
JSR_rep = 3;

%% ===== Modulos =====
sig = ti_signals();
noi = ti_noise();
itf = ti_interference();

%% ===== Tablas resultados =====
Tfull = table();
Trep  = table();

%% ===== Validacion NBI con f0 fijo =====
if (f0 - Bj/2) <= 0 || (f0 + Bj/2) >= (cfgBase.Fs/2)
    error("Bj=%.1f Hz no es valido para f0=%.1f Hz con Fs=%.1f. Ajusta Bj o f0.", Bj, f0, cfgBase.Fs);
end

%% ===== Bucle principal =====
for si = 1:numel(signals)
    sigType = signals{si};

    % Señal limpia x
    [~, xrow] = sig.generate(sigType, cfgBase);
    x = normalize_to_unit_power(remove_dc(xrow(:)));

    % PSD referencia para KL
    p_ref = psd_prob_from_welch(x, cfgBase.Fs, Wcfg, epsPSD);

    for ni = 1:numel(noises)
        noiseType = noises{ni};

        % Ruido (una vez por combinacion signal+noise)
        cfgN = cfgBase;
        cfgN.SNRdB  = SNRdB;
        cfgN.seed   = 1000 + 100*si + 10*ni;
        cfgN.useDSP = true;
        [~, wrow] = noi.generate(noiseType, cfgN);
        w = wrow(:);

        % ===== Baseline: SOLO RUIDO (sin interferencia) =====
        y0 = x + w;

        p0  = psd_prob_from_welch(y0, cfgBase.Fs, Wcfg, epsPSD);
        H0  = shannon_norm(p0);
        KL0 = kl_div(p_ref, p0);
        [MI0, ~, HYgX0, Eff0] = mi_and_derived(x, y0, MIcfg);

        % ===== Barrido JSR con 3 interferencias =====
        nJ  = numel(JSRdB_list);
        nIt = numel(itfList);

        Hobs = nan(nJ,nIt);
        dH   = nan(nJ,nIt);
        KL   = nan(nJ,nIt);
        MI   = nan(nJ,nIt);
        HYgX = nan(nJ,nIt);
        Eff  = nan(nJ,nIt);

        for jj = 1:nJ
            JSRdB = JSRdB_list(jj);

            for ii = 1:nIt
                itfType = itfList{ii};

                cfgJ = cfgBase;
                cfgJ.JSRdB  = JSRdB;
                cfgJ.seed   = 2000 + 100*si + 10*ni + 7*jj + ii;
                cfgJ.useDSP = true;

                switch itfType
                    case 'cw'
                        cfgJ.fj = f0;
                    case 'nbi'
                        cfgJ.fj = f0;
                        cfgJ.Bj = Bj;
                    case 'chirp'
                        cfgJ.f1 = max(10, f0 - chirpDelta);
                        cfgJ.f2 = min(cfgBase.Fs/2 - 10, f0 + chirpDelta);
                end

                [~, jrow] = itf.generate(itfType, cfgJ);
                j = jrow(:);

                y = x + w + j;

                % Espectrales
                p_obs = psd_prob_from_welch(y, cfgBase.Fs, Wcfg, epsPSD);
                Hobs(jj,ii) = shannon_norm(p_obs);
                dH(jj,ii)   = Hobs(jj,ii) - H0;
                KL(jj,ii)   = kl_div(p_ref, p_obs);

                % Temporales
                [MI(jj,ii), ~, HYgX(jj,ii), Eff(jj,ii)] = mi_and_derived(x, y, MIcfg);

                % ---- fila full (incluye baseline solo ruido) ----
                Tfull = [Tfull; table( ...
                    string(sigType), string(noiseType), string(itfType), ...
                    f0, SNRdB, JSRdB, ...
                    dH(jj,ii), MI(jj,ii), HYgX(jj,ii), Eff(jj,ii), KL(jj,ii), ...
                    H0, MI0, HYgX0, Eff0, KL0, ...
                    'VariableNames', { ...
                        'signal','noise','interference', ...
                        'f0_Hz','SNRdB','JSRdB', ...
                        'dH','MI','H_Y_given_X','Eff','KL', ...
                        'H_noiseonly','MI_noiseonly','H_Y_given_X_noiseonly','Eff_noiseonly','KL_noiseonly' ...
                    } ...
                )]; %#ok<AGROW>
            end
        end

        % ---- CSV resumen en JSR_rep (incluye baseline solo ruido) ----
        [~, idxRep] = min(abs(JSRdB_list - JSR_rep));
        for ii = 1:nIt
            dH_rep = dH(idxRep,ii);
            Trep = [Trep; table( ...
                string(sigType), string(noiseType), string(itfList{ii}), ...
                f0, SNRdB, JSRdB_list(idxRep), ...
                dH_rep, MI(idxRep,ii), HYgX(idxRep,ii), Eff(idxRep,ii), KL(idxRep,ii), ...
                H0, MI0, HYgX0, Eff0, KL0, ...
                'VariableNames', { ...
                    'signal','noise','interference', ...
                    'f0_Hz','SNRdB','JSRdB', ...
                    'dH','MI','H_Y_given_X','Eff','KL', ...
                    'H_noiseonly','MI_noiseonly','H_Y_given_X_noiseonly','Eff_noiseonly','KL_noiseonly' ...
                } ...
            )]; %#ok<AGROW>
        end

        % ---- Figura por (señal, ruido) ----
        fig = figure('Name', sprintf('%s + %s | f0=%.0f Hz, SNR=%.0f dB', sigType, noiseType, f0, SNRdB), ...
                     'Position',[100 100 1300 800]);

        ax1 = subplot(2,3,1);
        h1 = plot_multi(ax1, JSRdB_list, dH, itfList);
        grid(ax1,'on'); xlabel(ax1,'JSR (dB)');
        title(ax1, 'Variación Entropia  H(x+w+j) - H(x+w)', 'Interpreter','tex');
        yline(ax1, 0, '--');
        legend(ax1, h1, itfList, 'Location','best');

        ax2 = subplot(2,3,2);
        h2 = plot_multi(ax2, JSRdB_list, MI, itfList);
        grid(ax2,'on'); xlabel(ax2,'JSR (dB)');
        title(ax2, 'I(X;Y)', 'Interpreter','tex');
        legend(ax2, h2, itfList, 'Location','best');

        ax3 = subplot(2,3,3);
        h3 = plot_multi(ax3, JSRdB_list, HYgX, itfList);
        grid(ax3,'on'); xlabel(ax3,'JSR (dB)');
        title(ax3, 'H(Y|X) (ruido del canal)', 'Interpreter','tex');
        legend(ax3, h3, itfList, 'Location','best');

        ax4 = subplot(2,3,4);
        h4 = plot_multi(ax4, JSRdB_list, Eff, itfList);
        grid(ax4,'on'); xlabel(ax4,'JSR (dB)');
        title(ax4, 'Eff = I/H(Y)', 'Interpreter','tex');
        legend(ax4, h4, itfList, 'Location','best');

        ax5 = subplot(2,3,5);
        h5 = plot_multi(ax5, JSRdB_list, KL, itfList);
        grid(ax5,'on'); xlabel(ax5,'JSR (dB)');
        title(ax5, 'KL(p_x || p_y)', 'Interpreter','tex');
        legend(ax5, h5, itfList, 'Location','best');

        ax6 = subplot(2,3,6);
        axis(ax6,'off');
        text(ax6, 0, 0.85, sprintf('Señal=%s\nRuido=%s\nf0=%.1f Hz\nSNR=%.1f dB', ...
            sigType, noiseType, f0, SNRdB), 'FontSize', 12);

        text(ax6, 0, 0.35, sprintf('Señal limpia vs señal ruidosa:\nH=%.3f  MI=%.3f\nH(Y|X)=%.3f  Eff=%.3f\nKL=%.3f', ...
            H0, MI0, HYgX0, Eff0, KL0), 'FontSize', 10);

        drawnow;

        % Guarda figura
        figPath = fullfile(outDir, sprintf('fig_jsr_%s_%s.png', sigType, noiseType));
        exportgraphics(fig, figPath);
        close(fig);
    end
end

%% ===== Guardar CSVs =====
fullPath = fullfile(outDir, 'csv_full_jsr_sweep.csv');
repPath  = fullfile(outDir, 'csv_rep_jsr.csv');
writetable(Tfull, fullPath);
writetable(Trep,  repPath);

fprintf('Listo:\n  %s\n  %s\n', fullPath, repPath);
fprintf('Incluye baseline en columnas *_noiseonly y dH = H(with jam) - H(noise-only)\n');

end

%% ===================== helpers =====================

function h = plot_multi(ax, x, Y, labels)
axes(ax); %#ok<LAXES>
hold(ax, 'on');
h = gobjects(size(Y,2),1);
for i = 1:size(Y,2)
    h(i) = plot(ax, x, Y(:,i), 'LineWidth', 1.5, 'DisplayName', labels{i});
end
hold(ax, 'off');
end

function x = remove_dc(x)
x = x - mean(x);
end

function x = normalize_to_unit_power(x)
Px = mean(x.^2);
if Px <= 0, Px = 1; end
x = x / sqrt(Px);
end

function p = psd_prob_from_welch(x, Fs, Wcfg, epsPSD)
[Pxx, ~] = pwelch(x(:), Wcfg.winLen, Wcfg.overlap, Wcfg.nfft, Fs, 'onesided');
Pxx = max(Pxx(:), epsPSD);
p = Pxx / sum(Pxx);
end

function Hn = shannon_norm(p)
p = max(p(:), 1e-12);
p = p / sum(p);
M = numel(p);
H = -sum(p .* log2(p));
Hn = H / log2(M);
end

function D = kl_div(P, Q)
P = max(P(:), 1e-12); P = P / sum(P);
Q = max(Q(:), 1e-12); Q = Q / sum(Q);
D = sum(P .* log2(P ./ Q));
end

function [MI, HY, HYgX, Eff] = mi_and_derived(x, y, cfg)
x = x(:); y = y(:);

N = numel(x);
Neff = min(N, cfg.Neff_max);
idx = randperm(N, Neff);
x = x(idx); y = y(idx);

L = cfg.L; B = cfg.B;
edges = linspace(-L, L, B+1);

x = min(max(x, -L), L);
y = min(max(y, -L), L);

Nxy = histcounts2(x, y, edges, edges);
pXY = Nxy / max(1, sum(Nxy(:)));

pX = sum(pXY, 2);
pY = sum(pXY, 1).';

eps0 = cfg.eps;
pXY = max(pXY, eps0); pXY = pXY / sum(pXY(:));
pX  = max(pX,  eps0); pX  = pX  / sum(pX);
pY  = max(pY,  eps0); pY  = pY  / sum(pY);

HX  = -sum(pX .* log2(pX));
HY  = -sum(pY .* log2(pY));
HXY = -sum(pXY(:) .* log2(pXY(:)));

MI = HX + HY - HXY;
HYgX = HY - MI;

if HY > 0
    Eff = MI / HY;
else
    Eff = NaN;
end
end


