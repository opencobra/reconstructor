function [prefixMySQL,database] = initialiseMySqlCommand(inFile)
% initialiseMySqlCommand - Container-friendly VMH MySQL command prefix.
%
% This shadows the VMH toolbox helper from /matlab/scripts so MATLAB code can
% connect to a host/container MySQL server over TCP instead of assuming a Unix
% socket exists inside the MATLAB container.

if exist('inFile','var')
    infile = '--local-infile=1';
else
    infile = '';
end

database = getenvRequired('VMH_MYSQL_DATABASE');
host = getenvDefault('VMH_MYSQL_HOST', '127.0.0.1');
port = getenvDefault('VMH_MYSQL_PORT', '3306');
user = getenvDefault('VMH_MYSQL_USER', 'saleh');
password = getenv('VMH_MYSQL_PASSWORD');
extraArgs = getenv('VMH_MYSQL_EXTRA_ARGS');

parts = {};
if ~isempty(password)
    % MYSQL_PWD avoids the mysql client warning that can pollute captured output.
    parts{end + 1} = ['MYSQL_PWD=', shellQuote(password)];
end
parts{end + 1} = 'mysql';
if ~isempty(infile)
    parts{end + 1} = infile;
end
parts = [parts, {'--table', '--protocol=tcp', '-h', shellQuote(host), '-P', shellQuote(port), '-u', shellQuote(user)}];
if ~isempty(extraArgs)
    parts{end + 1} = extraArgs;
end
parts{end + 1} = shellQuote(database);

prefixMySQL = [strjoin(parts, ' '), ' -e'];
end

function value = getenvDefault(name, defaultValue)
value = getenv(name);
if isempty(value)
    value = defaultValue;
end
end

function value = getenvRequired(name)
value = getenv(name);
if isempty(value)
    error('Required environment variable %s is not set.', name);
end
end

function quoted = shellQuote(value)
value = char(value);
quoted = ['''', strrep(value, '''', '''"''"'''), ''''];
end
