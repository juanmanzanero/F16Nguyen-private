%
% Re-write the aero-dataset to add "numerical" breakpoints
% around beta = 0, so that the autodiff captures:
% 
%   * Zero slope at beta = 0 for the "CX", "CZ", "CM" coeffs.
%
%   * Symmetric slope at beta = 0, and value = 0 at beta = 0,
%     for the "CLL", "CY", "CM" coeffs.
%

outdir = 'aero_betasym';


% copy the original aero into a new aero folder,
% check that it is a new dir, or different to the
% original one
diroutdir   = dir(outdir);
diroriginal = dir('./aero');

assert(~exist(cell2mat(unique({diroutdir.folder})), 'dir') || ...
    ~strcmp(unique({diroutdir.folder}), unique({diroriginal.folder})));

system(['mkdir -p ', outdir]);
system(['cp -R ./aero/* ', outdir]);


% rewrite it with the new aos breaks
Dbeta_numerical_deg = 1e-4;
new_aos_breaks_deg  = sort([-Dbeta_numerical_deg; Dbeta_numerical_deg; ...
    read_dataset_csv([outdir, '/aos_breaks_deg.csv'])]);


% correct every table that depends on the aos
for d = diroutdir'
    if d.isdir
        continue
    end
    
    filename    = [d.folder, '/', d.name];
    [mat, info] = read_dataset_csv(filename);

    if contains(filename, 'aos')       
        if contains(filename, 'breaks')
            % breakpoints file
            assert(isequal(cell2mat(fieldnames(info)), 'title'))
            
            write_dataset_csv_breakpoints(filename, info.title, ...
                'aos_breaks_deg', new_aos_breaks_deg, '%.10f')

        else
            % table f(aos), add 2 new columns around beta = 0,
            % with the same values as that one (except for the
            % "CY", "CLL", and "CN" coefficients, which are not
            % symmetric w.r.t. beta normally)
            assert(strcmp(info.column_names, 'aos_breaks_deg'))

            beta0_column = info.columns == 0;
            assert(sum(beta0_column) == 1);
            pos_beta0_column = find(beta0_column);

            if contains(info.title{2}, 'CY') || ...
                    contains(info.title{2}, 'CN') || ...
                    contains(info.title{2}, 'CLL')
                % asymmetric coeffs: constant slope, and the
                % static ones are zero at beta = 0
                central_columns = interp1([info.columns(pos_beta0_column - 1) info.columns(pos_beta0_column + 1)], ...
                    [mat(:, pos_beta0_column - 1), mat(:, pos_beta0_column + 1)]', ...
                    [-Dbeta_numerical_deg, 0., Dbeta_numerical_deg])';
                                
                if contains(info.title{2}, 'CLL(') || ...
                        contains(info.title{2}, 'CY(') || ...
                        contains(info.title{2}, 'CN(')
                    
                    central_columns = central_columns - central_columns(:, 2);
                    
                end

                mat = [mat(:, 1:pos_beta0_column - 1), ...
                    central_columns, ...
                    mat(:, pos_beta0_column + 1:end)];

            elseif contains(info.title{2}, 'CX') || ...
                    contains(info.title{2}, 'CZ') || ...
                    contains(info.title{2}, 'CM')
                % symmetric coeffs: zero slope
                mat_beta0 = mat(:, pos_beta0_column);

                central_columns = [mat_beta0, mat_beta0, mat_beta0];
                
                mat = [mat(:, 1:pos_beta0_column - 1), ...
                    central_columns, ...
                    mat(:, pos_beta0_column + 1:end)];

            else
                error('rewrite_aero_with_beta_symmetry: wtf!')

            end

            write_dataset_csv_table(filename, info.title, ...
                mat, ...
                info.row_names, info.rows, ...
                info.column_names, new_aos_breaks_deg, '%.10f');


        end
        
    end
    
end


%
% Helper functions.
%

function [mat, info] = read_dataset_csv(filename)
%
% Reads a "*.csv" file from our datasets,
% as a matrix.
%

mat = readtable(filename, ...
    'ReadVariableNames', false, ...
    'FileType', 'text', ...
    'CommentStyle', '#');

mat = mat{:, :};


% the info are the title (first 2 commented lines) and the breaks
strfile    = fileread(filename);
title      = regexp(strfile, '(#.*?\n){2}', 'tokens', 'lineanchors', 'once');
rowinfo    = regexp(strfile, '\# ROWS (.*?)\:(.*?)$', 'tokens', 'lineanchors');
columninfo = regexp(strfile, '\# COLUMNS (.*?)\:(.*?)$', 'tokens', 'lineanchors');

if isempty(title) || isempty(title{:})
    error('rewrite_aero_with_beta_symmetry:read_dataset_csv: inconsistent table.')
else
    info = struct('title', {strtrim(strrep(strsplit(strtrim(title{:}), '\n'), '#', ''))});
end

if isempty(rowinfo) && isempty(columninfo)
    % a breakpoints file....
    return
    
else
    if ~isempty(rowinfo)
        num_rows = regexp(strfile, '\# NUM ROWS:(.*?)$', 'tokens', 'lineanchors');
        num_rows = str2double(num_rows{:}{:});
        
        info.row_names = rowinfo{1}{1};
        info.rows      = str2num(rowinfo{1}{2}); %#ok<ST2NM>
        
        if num_rows ~= size(mat, 1) || numel(info.rows) ~= num_rows || ~issorted(info.rows)
            error('rewrite_aero_with_beta_symmetry:read_dataset_csv: inconsistent table.')
        end
        
    else
        error('rewrite_aero_with_beta_symmetry:read_dataset_csv: inconsistent table.')
        
    end
    
    if ~isempty(columninfo)
        num_columns = regexp(strfile, '\# NUM COLUMNS:(.*?)$', 'tokens', 'lineanchors');
        num_columns = str2double(num_columns{:}{:});
        
        info.column_names = columninfo{1}{1};
        info.columns      = str2num(columninfo{1}{2}); %#ok<ST2NM>
        
        if num_columns ~= size(mat, 2) || numel(info.columns) ~= num_columns || ~issorted(info.columns)
            error('rewrite_aero_with_beta_symmetry:read_dataset_csv: inconsistent table.')
        end
        
    end
    
end

end


function write_dataset_csv_table(outfilename, title, ...
    mat, ...
    row_names, rows, ...
    column_names, columns, ...
    float_fmt)
%
% Writes a "*.csv" file with a 2-D table for our datasets.
%

if nargin < 8 || isempty(float_fmt)
    float_fmt = '%.5f';
end


if isrow(mat) || iscolumn(mat)
    mat         = mat(:);
    num_rows    = numel(mat);
    num_columns = 1;
    
    assert(~isempty(rows))
    assert(num_rows == numel(rows));
    
else
    num_rows    = size(mat, 1);
    num_columns = size(mat, 2);
    
    assert(~isempty(rows))
    assert(num_rows == numel(rows));
    assert(~isempty(columns))
    assert(num_columns == numel(columns));
    
end


fid = fopen(outfilename, 'w');
for line = 1:numel(title)
    fprintf(fid, '# %s\n', title{line});
end


fprintf(fid, '# NUM ROWS: %d\n', num_rows);
if ~isempty(row_names)
    fprintf(fid, ['# ROWS %s: ', float_fmt, ...
        repmat([',', float_fmt], 1, num_rows - 1),'\n'], row_names, rows(:)');
    
end

fprintf(fid, '# NUM COLUMNS: %d\n', num_columns);
if ~isempty(column_names)
    fprintf(fid, ['# COLUMNS %s: ', float_fmt, ...
        repmat([',', float_fmt], 1, num_columns - 1),'\n'], column_names, columns(:)');
    
end

fprintf(fid, [float_fmt, repmat([',', float_fmt], 1, num_columns - 1), '\n'], mat');
fclose(fid);

end


function write_dataset_csv_breakpoints(outfilename, title, ...
    name_breaks, breaks, ...
    float_fmt)
%
% Writes a "*.csv" file of breakpoints for our datasets.
%

if nargin < 5 || isempty(float_fmt)
    float_fmt = '%.5f';
end


assert(isrow(breaks) || iscolumn(breaks));

fid = fopen(outfilename, 'w');
for line = 1:numel(title)
    fprintf(fid, '# %s\n', title{line});
end
if ~isempty(name_breaks)
    fprintf(fid, '# %s\n', name_breaks);
end

fprintf(fid, '# NUM ROWS: %d\n', numel(breaks));
fprintf(fid, '# NUM COLUMNS: %d\n', 1);
fprintf(fid, [float_fmt, '\n'], breaks(:)');
fclose(fid);

end