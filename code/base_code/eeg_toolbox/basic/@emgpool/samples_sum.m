function s = samples_sum (pool)

% function s = samples_sum()
%
% Liefert die Anzahl an Daten die der Dateienpool
% pro Kanal enthaelt. Da alle Dateien pro Kanal
% gleich viele Daten enthalten, wird nur ein Wert
% zurueckgegeben.

s = sum(pool.samples);
