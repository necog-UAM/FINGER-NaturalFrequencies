function Figure_Object = plot_matches_redgreen(Matches)
% Plot the matches of the fingerprinting results
% "Matches" is a truth table, whose non-zero elements are the correlation of the vectors. Zero elements simply indicate a "not match". 
% Modified funtion of Javier Santonja.

% Define colors
color_match_diagonal = [0.2, 0.7, 0.3]; % Green for matches on the diagonal
color_match_off_diagonal = [0.9, 0.3, 0.2]; % Red for matches off the diagonal
color_no_match = [1 1 1]; % White for non-matches (optional)

% Create a matrix to store color values
color_matrix = zeros(size(Matches, 1), size(Matches, 2), 3);

% Assign colors based on the position and value
for i = 1:size(Matches, 1)
    for j = 1:size(Matches, 2)
        if Matches(i, j) ~= 0 % Match (non-zero value)
            if i == j
                color_matrix(i, j, :) = color_match_diagonal; % Diagonal match
            else
                color_matrix(i, j, :) = color_match_off_diagonal; % Off-diagonal match
            end
        else
            color_matrix(i, j, :) = color_no_match; % Non-match
        end
    end
end

% Plot the matrix using imagesc
Figure_Object = imagesc(color_matrix);

% Set axis properties
axis equal tight;
set(gca, 'XTick', [], 'YTick', []); % Hide tick marks

% Add a color legend for reference (Optional)
% Create a dummy colorbar to represent green and red colors
c = colorbar;
c.Ticks = [0 0.5 1];
c.TickLabels = {'No Match', 'Off-Diagonal Match', 'Diagonal Match'};
colormap([color_no_match; color_match_off_diagonal; color_match_diagonal]);

end
