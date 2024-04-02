function disp (pool)

% Overload of display method

disp('EMG-Pool: ');

for i = 1 : length(pool.file)
    disp(['   ' pool.file{i}]);
end
