function [D, r] = read (pool, chnnr, n, blocknr)

% function [D, r] = read(pool, chnnr, n, blocknr)
%
% Read channel data from pool of Neuro Explorer 
% emg-files.
%
% The data will be read in blocks of 'n'
% data values. The value 'blocknr' declares which 
% n-block to read.
%
%   chnnr - the channel whose data should be read
%       n - how many samples
% blocknr - the block that should be read

D = [];
r = -1;

if nargin ~= 4
    display('Need input arguments');
    return;
end

chnnr = fix(chnnr);

if chnnr < 1
    display('Channel number must be at least one');
    return;
end

chnpos = find(pool.channel == chnnr, 1) - 1;

if numel(chnpos) == 0
    display(['Channel ' num2str(chnnr) ' not stored in any file']);
    return;
end

n = fix(n);

if n < 1
    display('n must be at least 1');
    return;
end

blocknr = fix(blocknr);

if blocknr < 1
    display('Block number must be an integer greater one');
    return;
end

% Definitionen von Konstanten einlesen, die fuer das Auslesen
% der emg-Datei Header wichtig sind
nexdefines; 

% [D, r] = c_read(pool.file, pool.samples, chnnr, chnpos, length(pool.channel), n, blocknr, SIZE_EMG_FILEHEADER(pool.header), BYTES_PER_SAMPLE);

% Bestimme als Erstes die Position des ersten Samples
s_start = (blocknr - 1) * n + 1;

% Als naechstes den Index der dazugehoerigen Datei
k = 1;
s_sum_start = 0;
while k < length(pool.samples) + 1 && s_sum_start + pool.samples(k) < s_start
    s_sum_start = s_sum_start + pool.samples(k);
    k = k + 1;
end

if k == length(pool.samples) + 1
    %state = 1; % Ende erreicht
    return;
end

% Dann die Position des letzten zu lesenden Samples
s_stop = s_start + n;

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

% Lese jetzt die Daten aus
%s_stop
%s_start
%pool.samples
r = s_stop - s_start; % Anzahl Werte, die tatsaechlich ausgelesen werden
D = zeros(1, r);

i = 1;
while k < j

    % Datei oeffnen
    [emgfid, message] = fopen(pool.file{k}, 'r');

    if emgfid < 0
	display(['Fehler beim oeffnen von ' pool.file{k} ':' message]);
	return;
    end

    % Leseposition auf Daten setzen
    if 0 == pool.vpBit(k, chnpos + 1)
	% Die Daten liegen als float-Werte vor
	status = fseek(emgfid, SIZE_EMG_FILEHEADER(pool.header) + (s_start - s_sum_start - 1) * length(pool.channel) * BYTES_EMG_DATASAMPLE + chnpos * BYTES_EMG_DATASAMPLE, 'bof');

	if status < 0
	    display(ferror(emgfid));
	    fclose(emgfid);
	    return;
        end

	n_read = s_sum_start + pool.samples(k) - s_start + 1;

	D(i : i + n_read - 1) = fread(emgfid, n_read, SIZE_EMG_DATASAMPLE, (length(pool.channel) - 1) * BYTES_EMG_DATASAMPLE);
    else
	% Die Daten liegen als A/D-Werte vor
	status = fseek(emgfid, SIZE_EMG_FILEHEADER(pool.header) + (s_start - s_sum_start - 1) * length(pool.channel) * BYTES_EMG_ADSAMPLE + chnpos * BYTES_EMG_ADSAMPLE, 'bof');

	if status < 0
	    display(ferror(emgfid));
	    fclose(emgfid);
	    return;
        end

	n_read = s_sum_start + pool.samples(k) - s_start + 1;
	
	D(i : i + n_read - 1) = pool.vpBit(k, chnpos + 1) * fread(emgfid, n_read, SIZE_EMG_ADSAMPLE, (length(pool.channel) - 1) * BYTES_EMG_ADSAMPLE);
    end

    i = i + n_read;

    s_sum_start = s_sum_start + pool.samples(k);
    s_start = s_sum_start + 1; % Erstes Sample in neuer Datei
    k = k + 1;

    % Datei wieder schliessen
    fclose(emgfid);

end

% Datei j oeffnen
[emgfid, message] = fopen(pool.file{k}, 'r');

if emgfid < 0
    display(['Fehler beim oeffnen von ' pool.file{k} ':' message]);
    return;
end

% Leseposition auf Daten setzen 
if 0 == pool.vpBit(k, chnpos + 1)
    % Die Daten liegen als float-Werte vor
    status = fseek(emgfid, SIZE_EMG_FILEHEADER(pool.header) + (s_start - s_sum_start - 1) * length(pool.channel) * BYTES_EMG_DATASAMPLE + chnpos * BYTES_EMG_DATASAMPLE, 'bof');

    if status < 0
	display(ferror(emgfid));
	fclose(emgfid);
	return;
    end

    n_read = s_stop - s_start;

    D(i : i + n_read - 1) = fread(emgfid, n_read, SIZE_EMG_DATASAMPLE, (length(pool.channel) - 1) * BYTES_EMG_DATASAMPLE);
else
    % Die Daten liegen als A/D-Werte vor
    status = fseek(emgfid, SIZE_EMG_FILEHEADER(pool.header) + (s_start - s_sum_start - 1) * length(pool.channel) * BYTES_EMG_ADSAMPLE + chnpos * BYTES_EMG_ADSAMPLE, 'bof');

    if status < 0
	display(ferror(emgfid));
	fclose(emgfid);
	return;
    end

    n_read = s_stop - s_start;

    D(i : i + n_read - 1) = pool.vpBit(k, chnpos + 1) * fread(emgfid, n_read, SIZE_EMG_ADSAMPLE, (length(pool.channel) - 1) * BYTES_EMG_ADSAMPLE);
end

% Datei wieder schliessen
fclose(emgfid);

% Abschliessend die Daten auf die gewuenschte Einheit skalieren
D = pool.scale * D;
