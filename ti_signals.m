%[text] TI\_BASIC\_SIGNALS Generación y plot de señales básicas (tiempo).
%[text] Uso: m = ti\_basic\_signals(); cfg = struct('Fs',10000,'T',0.2,'f0',100,'A',1,'phi',0,'D',0.25,'Rc',1000,'seed',1); S = m.generate\_all(cfg); % devuelve struct con t y señales m.plot\_all(S); % pinta todas en una figura
function out = ti_signals()

out.generate      = @generate;
out.generate_all  = @generate_all;
out.plot_all      = @plot_all;
end

% ----------------- Generación individual -----------------

function [t, x] = generate(type, cfg)
% GENERATE  Genera una señal real discreta x[n] para un tipo dado.

% Defaults
if ~isfield(cfg,'A'),   cfg.A = 1; end
if ~isfield(cfg,'phi'), cfg.phi = 0; end

Fs  = cfg.Fs;
T   = cfg.T;
N   = round(Fs*T);
t   = (0:N-1)/Fs;

A   = cfg.A;
f0  = cfg.f0;
phi = cfg.phi;

type = lower(strtrim(type));

switch type
    case {'sin'}
        x = A*cos(2*pi*f0*t + phi);

    case {'sq'}
        x = sign(sin(2*pi*f0*t + phi));
        x(x==0) = 1;      % sgn(0)=+1
        x = A*x;

    case {'rect'}
        if ~isfield(cfg,'D'), error("Rectangular requiere cfg.D en (0,1)."); end
        D = cfg.D;
        if ~(D>0 && D<1), error("cfg.D debe estar en (0,1)."); end

        phase = mod(f0*t, 1);    % frac(f0*t)
        x = A*ones(size(t));
        x(phase >= D) = -A;

    case {'tri'}
        x = (2*A/pi)*asin(sin(2*pi*f0*t + phi));

    case {'saw'}
        frac = mod(f0*t, 1);
        x = A*(2*frac - 1);

    case {'dsss'}
    if ~isfield(cfg,'Rc'), error("DSSS requiere cfg.Rc (chips/seg)."); end
    if ~isfield(cfg,'seed'), cfg.seed = 1; end

    Rc = cfg.Rc;
    Nc_real = Fs/Rc;
    Nc = round(Nc_real);
    if abs(Nc_real - Nc) > 1e-12
        error("Fs/Rc debe ser entero. Ajusta Rc para que Nc=Fs/Rc sea entero.");
    end

    rng(cfg.seed);
    K = ceil(N/Nc);
    a = 2*randi([0 1], K, 1) - 1;   % chips ±1 (columna)
    c = repelem(a, Nc);             % columna
    c = c(1:N);                     % columna Nx1

    carrier = cos(2*pi*f0*t(:) + phi);  % columna Nx1
    x = A*(c .* carrier);               % columna Nx1

    x = x.';  % si quieres devolver fila para plotear (opcional)

    otherwise
        error("Tipo '%s' no reconocido.", type);
end

end

% ----------------- Generación de todas -----------------

function S = generate_all(cfg)
% GENERATE_ALL  Genera todas las señales básicas y las devuelve en un struct.

types = {'sin','sq','rect','tri','saw','dsss'};
names = {'sin','sq','rect','tri','saw','dsss'};

S = struct();
[t, xsin] = generate('sin', cfg);
S.t = t;
S.Fs = cfg.Fs;
S.f0 = cfg.f0;

S.(names{1}) = xsin;

for k = 2:numel(types)
    [~, x] = generate(types{k}, cfg);
    S.(names{k}) = x;
end

end

% ----------------- Plot conjunto -----------------

function fig = plot_all(S)
% PLOT_ALL  Pinta todas las señales en una figura (subplots).

t = S.t;
fields = {'sin','sq','rect','tri','saw','dsss'};
titles = {'Sinusoidal','Cuadrada','Rectangular','Triangular','Diente de sierra','DSSS/GNSS (simpl.)'};

% Para que se vea limpio, mostramos un número limitado de muestras:
nPeriods = 7;                        % <-- pon 2 o 3
Tplot = nPeriods / S.f0;             % segundos a mostrar
maxSamples = min(numel(t), floor(S.Fs * Tplot));
maxSamples = max(maxSamples, 10);  
fig = figure('Name','Señales básicas (tiempo)');
for i = 1:numel(fields)
    subplot(3,2,i);
    x = S.(fields{i});
    plot(t(1:maxSamples), x(1:maxSamples), 'LineWidth', 1); 
    grid on;
    title(sprintf('%s (f0=%.1f Hz)', titles{i}, S.f0));
    xlabel('t (s)');
    ylabel('x[n]');
        pad = 0.15;
    ymin = min(x(1:maxSamples));
    ymax = max(x(1:maxSamples));
    yr = ymax - ymin;
    if yr == 0, yr = 1; end
    ylim([ymin - pad*yr, ymax + pad*yr]);
end

end

%[appendix]{"version":"1.0"}
%---
