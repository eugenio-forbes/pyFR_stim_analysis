classdef emgpool	

     properties (SetAccess = 'private', GetAccess = 'private')
          file     = {}; % Dateinamen
	  samples  = []; % Anzahl Samples pro Kanal
	  channel  = []; % Gemeinsame Kanaele
	  sampling = []; % Samplingraten
	  vpBit    = []; % Matrix mit VoltPerBit-Eintraegen 
	                 % (i-te Zeile == i-te Datei, j-te Spalte == j-ter gespeicherter Kanal)
	  header   = 0;  % Headerversion
	  scale    = 1;  % Skalierungsfaktor fuer Einheit 'Volt'
     end

     properties (SetAccess = 'private', GetAccess = 'public')
	  valid    = 0;  % Fehler aufgetreten
     end

     properties (SetAccess = 'private', Dependent = true)
	  limit;         % Maximaler Samplewert
     end

     methods
          %
	  % Zugriff auf 'limit'
	  %
	  function l = get.limit(pool)
	     l = realmax('double'); % Maximaler Samplewert
	  end
	  %
          % Konstruktor definieren
	  %
          function pool = emgpool(filenames, si_unit)
	     
	     if nargin < 1
		 display('Liste mit Dateinamen erwartet');
		 return;
	     end

	     n = length(filenames);

	     if n == 0
		 display('Liste mit Dateinamen erwartet');
		 return;
	     end

	     if nargin == 2
                 switch( si_unit )
		    case 'Volt'
		        pool.scale = 1;
		    case 'Millivolt'
		        pool.scale = 1e3;
		    case 'Mikrovolt'
		        pool.scale = 1e6;
		    otherwise
		        display('Unknown SI unit.');
	         end
	     end

	     % Laenge der Attributvektoren festlegen
	     pool.samples  = zeros(1, n);
	     pool.sampling = zeros(1, n);

	     % Definitionen von Konstanten einlesen, die fuer das Auslesen
	     % der emg-Datei Header wichtig sind
	     nexdefines; 

	     % Jede angegebene Datei oeffnen und pruefen, ob in allen
	     % die gleichen Kanaele gespeichert sind
	     for i = 1 : n

		 % emg Datei oeffnen
		 [emgfid, message] = fopen(filenames{i}, 'r');

		 if emgfid < 0
		     display(['Fehler beim oeffnen von ' filenames{i} ': ' message]);
		     return;
		 end

		 % An das Ende der emg-Datei springen und Position (= Dateigroesse) auslesen
		 status = fseek(emgfid, 0, 'eof');

		 if status < 0
		     display(ferror(emgfid));
		     fclose(emgfid);
		     return;
		 end

		 position = ftell(emgfid);

		 if position < 0
		     display(ferror(emgfid));
		     fclose(emgfid);
		     return;
		 end

		 if position < POSITION_HEADER_VERSION_INFO + 4
		     display(['File' filenames{i} ' does not contain any data or is incomplete']);
		     fclose(emgfid);
		     return;
		 end

		 % Version des Header auslesen
		 status = fseek(emgfid, POSITION_HEADER_VERSION_INFO, 'bof');

		 if status < 0
		     display(ferror(emgfid));
		     fclose(fid);
		     return;
		 end

		 tmp_hv = fread(emgfid, 1, SIZE_HEADER_VERSION_INFO) + 1;

		 if pool.header == 0
		     pool.header = tmp_hv;
		 else
		     if pool.header ~= tmp_hv
			 display(['Datei ' filenames{i} ' hat eine andere Headergroesse als vorhergehende Dateien']);
			 fclose(emgfid);
			 return;
		     end
		 end

		 % Nochmal pruefen, ob Daten enthalten sind
		 if position < SIZE_EMG_FILEHEADER(tmp_hv)
		     display(['File' filenames{i} ' does not contain any data or is incomplete']);
		     fclose(emgfid);
		     return;
		 end

		 % Anzahl gespeicherter Kanaele auslesen
		 status = fseek(emgfid, POSITION_EMG_NUMCHANNEL, 'bof');

		 if status < 0
		     display(ferror(emgfid));
		     fclose(emgfid);
		     return;
		 end

		 nCh = fread(emgfid, 1, SIZE_EMG_NUMCHANNEL);

		 % Lese aus, welche Kanaele gespeichert sind sowie den VoltPerBit-Eintrag
		 % jedes gespeicherten Kanals
		 tmp_chn = [];
		 tmp_vpBit = [];
		 active  = 0;
		 k = 0;
		 while k < NUM_EMG_CHANNELS(tmp_hv)
		     status = fseek(emgfid, POSITION_EMG_CHANNELINFO + k * SIZE_CHANNEL_STRUCT + POSITION_CHANNELFLAG, 'bof');

		     if status < 0
			 display(ferror(emgfid));
			 fclose(emgfid);
			 return;
		     end

		     active = fread(emgfid, 1, SIZE_CHANNELFLAG);

		     if 1 == bitand(active, 1)
			 tmp_chn(end + 1) = k + 1;

			 status = fseek(emgfid, POSITION_EMG_CHANNELINFO + k * SIZE_CHANNEL_STRUCT + POSITION_VOLTPERBIT, 'bof');

			 if status < 0
			     display(ferror(emgfid));
			     fclose(emgfid);
			     return;
		         end

		         tmp_vpBit(end + 1) = fread(emgfid, 1, SIZE_VOLTPERBIT);
		     end

		     k = k + 1;
		 end

		 if length(tmp_chn) == 0
		     display(['No data stored in file']);
		     fclose(emgfid);
		     return;
		 end

		 if prod(size(pool.channel)) == 0
		     pool.channel = tmp_chn;
		 else
		     if (length(pool.channel) ~= length(tmp_chn) || (length(pool.channel) == length(tmp_chn) && any(pool.channel - tmp_chn)))
			 display(['Datei ' filenames{i} ' enthaelt andere Kanaele als vorhergehende Dateien']);
			 fclose(emgfid);
			 return;
		     end
		 end

		 % Zusammen mit ermittelter Dateigroesse die Anzahl an Samples pro Kanal bestimmen
		 ze = sum(0 == tmp_vpBit);

		 if 0 == ze
		     % Alle VoltPerBit-Eintraege ungleich 0, d.h. es liegt eine Datei mit A/D-Werten vor
		     pool.samples(i) = ((position - SIZE_EMG_FILEHEADER(tmp_hv)) / BYTES_EMG_ADSAMPLE ) / nCh;
		 elseif ze == length(tmp_vpBit)
		     % Alle VoltPerBit-Eintraege gleich 0, d.h. es liegt eine Datei mit float-Werten vor
		     pool.samples(i) = ((position - SIZE_EMG_FILEHEADER(tmp_hv)) / BYTES_EMG_DATASAMPLE ) / nCh;
	         else
		     % Einige VoltPerBit-Eintraege sind gleich 0, einige ungleich 0, d.h. mit der Datei stimmt was nicht
		     display(['Datei ' filenames{i} ' enthaelt merkwuerdige VoltPerBit-Eintraege']);
		     fclose(emgfid);
		     return;
	         end

		 % Gesammelte Eintraege anhaengen
		 pool.vpBit = cat(2, pool.vpBit, tmp_vpBit);

		 % Read sampling rate
		 status = fseek(emgfid, POSITION_EMG_SAMPLING, 'bof');

		 if status < 0
		     display(ferror(fid));
		     fclose(fid);
		     return;
		 end

		 pool.sampling(i) = 20000;%fread(emgfid, 1, SIZE_EMG_SAMPLING);

		 if pool.sampling(i) == 0
		     display(['Zero sampling rate in file' filenames{i}]);
		     fclose(emgfid);
		     return;
		 end

		 % emg Datei wieder schliessen
		 fclose(emgfid);

		 %
	     end

	     % Matrix mit VoltPerBit-Eintraegen aus Vektor erzeugen
	     pool.vpBit = reshape(pool.vpBit, length(pool.channel), n)';

	     % Dateinamen merken
	     pool.file = filenames;

	     % Konstruktion gelungen
	     pool.valid = 1;

          end
     end
end
