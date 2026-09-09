%% ========================================================================
%  CHECK_PACKAGES  --  is every file where MATLAB expects it?
%
%      >> cd C:\path\to\proj1
%      >> check_packages
%
%  Catches the two things that keep going wrong when files are downloaded
%  one at a time: a lost or doubled .m extension, and a stale path cache
%  after adding a file to an existing +package folder.
% ========================================================================
clc;
clearvars

fprintf('========================================================\n');
fprintf(' package integrity check\n');
fprintf('========================================================\n\n');
fprintf('Folder: %s\n\n', pwd);

expected = struct( ...
    'windlib',    {{'load','registry','resample','feasibility'}}, ...
    'predictors', {{'registry','available','cached','num_segments','stitch','rbf_core', ...
                    'rbf','rbf_arima','arima','svr','gp','esn','gru','lstm','bilstm','tcn'}}, ...
    'ddp',        {{'params','refs','solve','forward','backward','rk4', ...
                    'second_order','compute_Ta_dot','resample_segment'}}, ...
    'report',     {{'paths','build_data','update_manifest','export_prediction', ...
                    'export_control','segment_figures','stats_figures','animation','write_txt'}} );

pkgs = fieldnames(expected);
problems = 0;

rehash;   % force a re-scan before we look

for i = 1:numel(pkgs)
    pkg = pkgs{i};
    dirname = ['+' pkg];
    fprintf('--- %s ---\n', dirname);

    if exist(dirname, 'dir') ~= 7
        fprintf(2, '  MISSING FOLDER: %s\n', fullfile(pwd, dirname));
        problems = problems + 1;
        fprintf('\n');
        continue;
    end

    files = expected.(pkg);
    listing = dir(dirname);
    onDisk = {listing(~[listing.isdir]).name};

    missing = {}; unresolved = {};
    for j = 1:numel(files)
        fname = [files{j} '.m'];

        if ~any(strcmp(onDisk, fname))
            % look for a mangled variant before giving up
            hint = '';
            for k = 1:numel(onDisk)
                [~, base, ext] = fileparts(onDisk{k});
                if strcmpi(base, files{j}) || strcmpi(onDisk{k}, [files{j} '.m.txt'])
                    hint = sprintf('  (found "%s" -- rename it to "%s")', onDisk{k}, fname);
                    break;
                end
            end
            missing{end+1} = [files{j} hint]; %#ok<SAGROW>
            continue;
        end

        % present on disk -- can MATLAB actually resolve it?
        if isempty(which([pkg '.' files{j}]))
            unresolved{end+1} = files{j}; %#ok<SAGROW>
        end
    end

    if isempty(missing) && isempty(unresolved)
        fprintf('  OK  all %d files present and resolvable\n', numel(files));
    end
    for j = 1:numel(missing)
        fprintf(2, '  MISSING     %s\n', missing{j});
        problems = problems + 1;
    end
    for j = 1:numel(unresolved)
        fprintf(2, '  UNRESOLVED  %s.m is on disk but MATLAB cannot see it\n', unresolved{j});
        fprintf(2, '              try:  rehash toolboxcache\n');
        problems = problems + 1;
    end

    % anything unexpected sitting in the package folder
    extra = setdiff(onDisk, strcat(files, '.m'));
    extra = extra(~strcmp(extra, 'data'));
    if ~isempty(extra)
        fprintf('  note: also present -> %s\n', strjoin(extra, ', '));
    end
    fprintf('\n');
end

%% ---- data files ----
fprintf('--- +windlib/data ---\n');
dd = fullfile('+windlib','data');
if exist(dd,'dir') ~= 7
    fprintf(2, '  MISSING FOLDER: %s\n', fullfile(pwd, dd));
    if exist('data','dir') == 7
        fprintf(2, '  There is a data folder at %s -- move it INSIDE +windlib.\n', fullfile(pwd,'data'));
    end
    problems = problems + 1;
else
    csvs = {'wind01','wind02','wind03','wind06','wind07'};
    miss = 0;
    for i = 1:numel(csvs)
        if exist(fullfile(dd,[csvs{i} '.csv']),'file') ~= 2
            fprintf(2, '  MISSING %s.csv\n', csvs{i});
            miss = miss + 1;
        end
    end
    if miss == 0
        fprintf('  OK  all 5 csv files present\n');
    end
    problems = problems + miss;
end
fprintf('\n');

%% ---- root scripts ----
fprintf('--- root scripts ---\n');
scripts = {'run_prediction','run_control','report.html'};
for i = 1:numel(scripts)
    s = scripts{i};
    if endsWith(s, '.html')
        ok = exist(s,'file') == 2;
    else
        ok = exist([s '.m'],'file') == 2;
    end
    if ok
        fprintf('  OK      %s\n', s);
    else
        fprintf(2, '  MISSING %s\n', s);
        problems = problems + 1;
    end
end
fprintf('\n');

%% ---- shadowing check ----
fprintf('--- variable shadowing ---\n');
shadow = 0;
for i = 1:numel(pkgs)
    if evalin('base', sprintf('exist(''%s'',''var'')', pkgs{i}))
        fprintf(2, '  A VARIABLE named "%s" is shadowing the package. Run: clear %s\n', ...
            pkgs{i}, pkgs{i});
        shadow = shadow + 1;
    end
end
if shadow == 0
    fprintf('  OK  no package name is shadowed by a variable\n');
end
problems = problems + shadow;
fprintf('\n');

%% ---- verdict ----
fprintf('========================================================\n');
if problems == 0
    fprintf(' ALL PACKAGES INTACT\n');
else
    fprintf(2, ' %d problem(s) found -- see above\n', problems);
end
fprintf('========================================================\n');
