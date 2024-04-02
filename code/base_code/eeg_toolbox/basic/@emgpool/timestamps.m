function [T, r] = timestamps (pool, n, blocksize, seqnr)

% function T = timestamps (n, blocksize, seqnr)
%
% Liefert fuer die Sequenz 'seqnr' Zeitstempel fuer 'n' aufeinanderfolgende 
% Datenbloecke mit je 'blocksize' Werten. Hierzu werden die Daten
% in Sequenzen der Laenge 'n * blocksize' unterteilt.
%
% T -- Vektor mit Zeitstempeln
% r -- Anzahl Zeitstempel, r <= n, oder
%      -1 bei Fehler oder r = 0, wenn
%      das Ende der Daten erreicht ist

T = [];
r = -1;

if nargin ~= 4
    display('Need input arguments');
    return;
end

n = fix(n);

if n < 1
    display('n must be at least 1');
    return;
end

blocksize = fix(blocksize);

if blocksize < 1
    display('Block size must be an integer greater one');
    return;
end

seqnr = fix(seqnr);

if seqnr < 1
    display('Sequence number must be at least one');
    return;
end

% Rufe C-Funktion auf
[T, r] = c_timestamps(pool.samples, pool.sampling, n, blocksize, seqnr);

% Skaliere auf µs
T = 1e6 * T;



%n_sample = length(pool.samples) + 1;
%
%s = (seqnr - 1) * n * blocksize + 1;
%
%k = 1; % In C dann von 0 beginnen
%s_sum = 0;
%t_off = 0;
%while k < n_sample && s_sum + pool.samples(k) < s 
%    t_off = t_off + pool.samples(k) / pool.sampling(k);
%    s_sum = s_sum + pool.samples(k);
%    k = k + 1;
%end
%
%if k == n_sample
%    r = 0; % Ende erreicht
%    return;
%end
%
%t_off = t_off + (s - s_sum - 1) / pool.sampling(k);
%
%% Koennen alle vorgesehenen Zeitstempel bestimmt werden?
%r = n;
%
%s_new = s + (n - 1) * blocksize;
%j = k;
%s_sum_new = s_sum;
%while j < n_sample && s_sum_new + pool.samples(j) < s_new 
%    s_sum_new = s_sum_new + pool.samples(j);
%    j = j + 1;
%end
%
%if j == n_sample
%    r = ceil((s_sum_new - s + 1) / blocksize);
%end
%
%T = zeros(1, r);
%T(1) = t_off;   % In C dann Index 0
%
%s_sum_old = s_sum;
%s_old = s;
%
%for i = 2 : r
%    s_new = s_old + blocksize; % Weiterruecken und pruefen, ob Dateigrenze ueberschritten
%
%    j = k;
%    s_sum_new = s_sum_old;
%    while j < n_sample && s_sum_new + pool.samples(j) < s_new 
%	s_sum_new = s_sum_new + pool.samples(j);
%	j = j + 1;
%    end
%
%    T(i) = T(i-1)
%    while k < j
%	T(i) = T(i) + (s_sum_old + pool.samples(k) - s_old + 1) / pool.sampling(k);
%	s_sum_old = s_sum_old + pool.samples(k);
%	s_old = s_sum_old + 1; % Erstes Sample in neuer Datei
%	k = k + 1;
%    end
%    T(i) = T(i) + (s_new - s_old) / pool.sampling(k);
%
%    s_old = s_new;
%end
