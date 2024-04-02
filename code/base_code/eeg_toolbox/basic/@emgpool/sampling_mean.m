function sfm = sampling_mean (pool, blocksize, blocknr)

% function sfm = sampling_mean (blocksize, blocknr)
%
% Liefert die mittlere Samplingrate fuer die Datenwerte
% des Blocks 'blocknr'. Die Daten werden hierzu in 
% Bloecke mit je 'blocksize' Werte unterteilt.

sfm = -1;

if nargin ~= 3
    display('Need input arguments');
    return;
end

blocksize = fix(blocksize);

if blocksize < 1
    display('Block size must be at least one');
    return;
end

blocknr = fix(blocknr);

if blocknr < 1
    display('Block number must be an integer greater one');
    return;
end

% Bestimme als erstes die Position des ersten Samples
s_start = (blocknr - 1) * blocksize + 1;

% Als naechstes den Index der dazugehoerigen Datei
k = 1;
s_sum_start = 0;
while k < length(pool.samples) + 1 && s_sum_start + pool.samples(k) < s_start
    s_sum_start = s_sum_start + pool.samples(k);
    k = k + 1;
end

if k == length(pool.samples) + 1
    sfm = -1; % Ende erreicht
    return;
end

% Dann die Position des letzten Samples
s_stop = s_start + blocksize;

% Und den Index der dazugehoerigen Datei
j = k;
s_sum_stop = s_sum_start;
while j < length(pool.samples) + 1 && s_sum_stop + pool.samples(j) < s_stop
    s_sum_stop = s_sum_stop + pool.samples(j);
    j = j + 1;
end

% Ende erreicht
if j == length(pool.samples) + 1
    j = length(pool.samples);
    s_stop = s_sum_stop + 1;
end

% Errechne jetzt die mittlere Samplingrate
n_sum = 0; % Summe der Samples
t_sum = 0; % Summe der verstrichenen Zeit

while k < j

    n = (s_sum_start + pool.samples(k)) - s_start + 1;

    n_sum = n_sum + n;
    t_sum = t_sum + n / pool.sampling(k);

    s_sum_start = s_sum_start + pool.samples(k);
    s_start     = s_sum_start + 1; % Erstes Sample in neuer Datei
    k = k + 1;

end

% Fuer j
n = s_stop - s_start;

n_sum = n_sum + n;
t_sum = t_sum + n / pool.sampling(j);

sfm = n_sum / t_sum;
