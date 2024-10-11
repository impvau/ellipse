
digits(200);

% List of specific a values (b is constant 1)
specific_a_values = [105/100, 115/100, 125/100, 135/100, 145/100, 155/100, 165/100, ...
                     175/100, 185/100, 195/100, 40, 50, 60, 70, 80, 90, 100, ...
                     500, 1000, 10000, 1000000, 1000000000];

% Generate a range of a values from 2 to 30
range_a_values = 2:30;

% Combine specific and range a values
a_values = [specific_a_values, range_a_values];

b = sym(1);  % Use symbolic precision for b

% Sort the a values in increasing order
a_values = sort(a_values);

% Open file to write results
fileID = fopen('elliptical_perimeter_results.txt', 'w');

% Compute and write perimeter for each value of a
for i = 1:length(a_values)
    a = sym(a_values(i)); % Use symbolic precision for a
    k = sqrt(1 - (b^2 / a^2)); % Calculate the modulus k
    elliptical_perimeter = vpa(4 * ellipticE(k^2) * a); % Use elliptic integral of the second kind with high precision
    
    % Write the result to the file with perimeter first, a second, separated by a comma
    fprintf(fileID, '%s, %s\n', char(elliptical_perimeter), char(a));
end

fclose(fileID);

disp('Elliptical perimeter values have been computed and saved to elliptical_perimeter_results.txt with 200 decimal places');
