function [RunCorrT, W, A, icasig, BestRun] = icatb_bestRunSelection(icasigR, data)
%% Select best run based on MST.
% Inputs:
% 1. icasigR - ICA signals in a cell array of length equal to the number of ICA runs
% 2. data - Data of dimensions components x voxels x subjects
%
% Outputs:
% RunCorrT - Correlation metric
% W - Unmixing matrix of best run
% A - Mixing matrix of best run
% icasig - Spatial sources of best run
% BestRun - Best run
%
% Program by Wei Du. Please contact me at weidu1@umbc.edu
% MLSP-Lab http://mlsp.umbc.edu


%% Initialize vars
numOfRuns = length(icasigR);
numOfIC = size(icasigR{1}, 1);
numOfVoxels = size(icasigR{1}, 2);
numOfSub = size(icasigR{1}, 3);

%% When using IVA stack subjects along voxel dimension
if (numOfSub > 1)
    for nR = 1:numOfRuns
        icasigR{nR} = reshape(icasigR{nR}, numOfIC, numOfVoxels*numOfSub);
    end
end

[Tmaps, IcSortMag] = TmapsMultiRuns(icasigR);

RunCorrT = zeros(numOfIC, numOfRuns);
for k = 1:numOfRuns
    RunCorrT(:, k) = diag(abs(icatb_corr(IcSortMag{k}', Tmaps')));
end

[~, BestRun] = max(sum(RunCorrT));

icasig = icasigR{BestRun};

if (numOfSub > 1)
    icasig = reshape(icasig, numOfIC, numOfVoxels, numOfSub);
    A = zeros(numOfIC, numOfIC, numOfSub);
    W = A;
    for nSub = 1:numOfSub
        A(:, :, nSub) = squeeze(data(:, :, nSub))*pinv(squeeze(icasig(:, :, nSub)));
        W(:, :, nSub) = pinv(A(:, :, nSub));
    end
else
    % Compute mixing and un-mixing matrices of ICA algorithm
    A = data*pinv(icasig);
    W = pinv(A);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Program by Wei Du. Please contact me at weidu1@umbc.edu
% MLSP-Lab http://mlsp.umbc.edu

function [Tmaps, icSort, IndexAssign] = TmapsMultiRuns(icasigR)
% % Components from each are clustered
[IndexAssign, icSort] = ComponentAlignment(icasigR);

NumRun = length(icasigR);
NumComp = size(icasigR{1}, 1);
voxels = size(icasigR{1}, 2);

for k = 1:NumRun
    icSort{k} = icatb_zscore(icSort{k}')';
end

Tmaps = zeros(NumComp, voxels);
for i = 1:NumComp
    icCrun = zeros(NumRun, voxels);
    for j = 1:NumRun
        icSrun = icSort{j};
        v = icSrun(i,:);
        v = icatb_recenter_image(v);
        skew = icatb_skewness(v);
        if (sign(skew) == -1)
            v = v*-1;
        end
        icCrun(j, :) = v;
    end
    [pval, tval] = icatb_ttest(icCrun);
    Tmaps(i, :) = tval;
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Program by Wei Du. Please contact me at weidu1@umbc.edu
% MLSP-Lab http://mlsp.umbc.edu

function [IndexAssign, icasigR] = ComponentAlignment(icasigR)

NumRun = length(icasigR);
NumComp = size(icasigR{1}, 1);

C = nchoosek(1:NumRun,2);
cost = zeros(1, size(C, 1));
assignment = cell(1, size(C, 1));
for i=1:size(C, 1)
    corrM = abs(icatb_corr(icasigR{1,C(i,1)}', icasigR{1,C(i,2)}'));
    corrM(corrM > 1) = 1;
    [assignment{1,i}, cost(1,i)] = assignmentoptimal(1 - corrM);
end

MatrixCost = squareform( cost);

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % MST
try
    [ST, pred] = graphminspantree(sparse(MatrixCost));
catch
    ST = icatb_mst(sparse(MatrixCost));
end
%view(biograph(ST,[],'ShowArrows','off','ShowWeights','on'))
[Treei,Treej] = find(ST);

Table_tmp = maketbl([Treei;Treej]);
[vijMax,tijmax] = max(Table_tmp(:,2));
TreeStart = Table_tmp(tijmax,1);


disp(['The central node is ',num2str(TreeStart)]);

IndexAssign = zeros(NumComp,NumRun);
IndexAssign(:,TreeStart) = (1:NumComp)';
tmp = squareform(1:size(C,1));

% % assign components according to the centrol node in MST
%
j=1;
while j <= NumRun
    if j~=TreeStart
        IndexAssign(:,j) = AssignPair(TreeStart, j, IndexAssign(:, TreeStart), tmp, assignment);
    end
    j=j+1;
end

for k = 1:NumRun
    tmpInd = IndexAssign(:, k);
    icasigR{k} = icasigR{k}(tmpInd, :);
end


function Resultij = AssignPair(indi,indj,Resulti,tmp,assignment)
if indj>indi
    Resultj=assignment{1,tmp(indi,indj)};
else
    tmp_assign = assignment{1,tmp(indi,indj)};
    [~,assignInd] = sort(tmp_assign);
    Resultj = assignInd;
end
Resultij = Resultj(Resulti);


function [assignment, cost] = assignmentoptimal(distMatrix)
%ASSIGNMENTOPTIMAL    Compute optimal assignment by Munkres algorithm
%		ASSIGNMENTOPTIMAL(DISTMATRIX) computes the optimal assignment (minimum
%		overall costs) for the given rectangular distance or cost matrix, for
%		example the assignment of tracks (in rows) to observations (in
%		columns). The result is a column vector containing the assigned column
%		number in each row (or 0 if no assignment could be done).
%
%		[ASSIGNMENT, COST] = ASSIGNMENTOPTIMAL(DISTMATRIX) returns the
%		assignment vector and the overall cost.
%
%		The distance matrix may contain infinite values (forbidden
%		assignments). Internally, the infinite values are set to a very large
%		finite number, so that the Munkres algorithm itself works on
%		finite-number matrices. Before returning the assignment, all
%		assignments with infinite distance are deleted (i.e. set to zero).
%
%		A description of Munkres algorithm (also called Hungarian algorithm)
%		can easily be found on the web.
%
%		<a href="assignment.html">assignment.html</a>  <a href="http://www.mathworks.com/matlabcentral/fileexchange/6543">File Exchange</a>  <a href="https://www.paypal.com/cgi-bin/webscr?cmd=_s-xclick&hosted_button_id=EVW2A4G2HBVAU">Donate via PayPal</a>
%
%		Markus Buehren
%		Last modified 05.07.2011

% save original distMatrix for cost computation
originalDistMatrix = distMatrix;

% check for negative elements
if any(distMatrix(:) < 0)
    error('All matrix elements have to be non-negative.');
end

% get matrix dimensions
[nOfRows, nOfColumns] = size(distMatrix);

% check for infinite values
finiteIndex   = isfinite(distMatrix);
infiniteIndex = find(~finiteIndex);
if ~isempty(infiniteIndex)
    % set infinite values to large finite value
    maxFiniteValue = max(max(distMatrix(finiteIndex)));
    if maxFiniteValue > 0
        infValue = abs(10 * maxFiniteValue * nOfRows * nOfColumns);
    else
        infValue = 10;
    end
    if isempty(infValue)
        % all elements are infinite
        assignment = zeros(nOfRows, 1);
        cost       = 0;
        return
    end
    distMatrix(infiniteIndex) = infValue;
end

% memory allocation
coveredColumns = zeros(1,       nOfColumns);
coveredRows    = zeros(nOfRows, 1);
starMatrix     = zeros(nOfRows, nOfColumns);
primeMatrix    = zeros(nOfRows, nOfColumns);

% preliminary steps
if nOfRows <= nOfColumns
    minDim = nOfRows;
    
    % find the smallest element of each row
    minVector = min(distMatrix, [], 2);
    
    % subtract the smallest element of each row from the row
    distMatrix = distMatrix - repmat(minVector, 1, nOfColumns);
    
    % Steps 1 and 2
    for row = 1:nOfRows
        for col = find(distMatrix(row,:)==0)
            if ~coveredColumns(col)%~any(starMatrix(:,col))
                starMatrix(row, col) = 1;
                coveredColumns(col)  = 1;
                break
            end
        end
    end
    
else % nOfRows > nOfColumns
    minDim = nOfColumns;
    
    % find the smallest element of each column
    minVector = min(distMatrix);
    
    % subtract the smallest element of each column from the column
    distMatrix = distMatrix - repmat(minVector, nOfRows, 1);
    
    % Steps 1 and 2
    for col = 1:nOfColumns
        for row = find(distMatrix(:,col)==0)'
            if ~coveredRows(row)
                starMatrix(row, col) = 1;
                coveredColumns(col)  = 1;
                coveredRows(row)     = 1;
                break
            end
        end
    end
    coveredRows(:) = 0; % was used auxiliary above
end

if sum(coveredColumns) == minDim
    % algorithm finished
    assignment = buildassignmentvector__(starMatrix);
else
    % move to step 3
    [assignment, distMatrix, starMatrix, primeMatrix, coveredColumns, coveredRows] = step3__(distMatrix, starMatrix, primeMatrix, coveredColumns, coveredRows, minDim); %#ok
end

% compute cost and remove invalid assignments
[assignment, cost] = computeassignmentcost__(assignment, originalDistMatrix, nOfRows);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function assignment = buildassignmentvector__(starMatrix)

[maxValue, assignment] = max(starMatrix, [], 2);
assignment(maxValue == 0) = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [assignment, cost] = computeassignmentcost__(assignment, distMatrix, nOfRows)

rowIndex   = find(assignment);
costVector = distMatrix(rowIndex + nOfRows * (assignment(rowIndex)-1));
finiteIndex = isfinite(costVector);
cost = sum(costVector(finiteIndex));
assignment(rowIndex(~finiteIndex)) = 0;

% Step 2: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [assignment, distMatrix, starMatrix, primeMatrix, coveredColumns, coveredRows] = step2__(distMatrix, starMatrix, primeMatrix, coveredColumns, coveredRows, minDim)

% cover every column containing a starred zero
maxValue = max(starMatrix);
coveredColumns(maxValue == 1) = 1;

if sum(coveredColumns) == minDim
    % algorithm finished
    assignment = buildassignmentvector__(starMatrix);
else
    % move to step 3
    [assignment, distMatrix, starMatrix, primeMatrix, coveredColumns, coveredRows] = step3__(distMatrix, starMatrix, primeMatrix, coveredColumns, coveredRows, minDim);
end

% Step 3: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [assignment, distMatrix, starMatrix, primeMatrix, coveredColumns, coveredRows] = step3__(distMatrix, starMatrix, primeMatrix, coveredColumns, coveredRows, minDim)

zerosFound = 1;
while zerosFound
    
    zerosFound = 0;
    for col = find(~coveredColumns)
        for row = find(~coveredRows')
            if distMatrix(row,col) == 0
                
                primeMatrix(row, col) = 1;
                starCol = find(starMatrix(row,:));
                if isempty(starCol)
                    % move to step 4
                    [assignment, distMatrix, starMatrix, primeMatrix, coveredColumns, coveredRows] = step4__(distMatrix, starMatrix, primeMatrix, coveredColumns, coveredRows, row, col, minDim);
                    return
                else
                    coveredRows(row)        = 1;
                    coveredColumns(starCol) = 0;
                    zerosFound              = 1;
                    break % go on in next column
                end
            end
        end
    end
end

% move to step 5
[assignment, distMatrix, starMatrix, primeMatrix, coveredColumns, coveredRows] = step5__(distMatrix, starMatrix, primeMatrix, coveredColumns, coveredRows, minDim);

% Step 4: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [assignment, distMatrix, starMatrix, primeMatrix, coveredColumns, coveredRows] = step4__(distMatrix, starMatrix, primeMatrix, coveredColumns, coveredRows, row, col, minDim)

newStarMatrix          = starMatrix;
newStarMatrix(row,col) = 1;

starCol = col;
starRow = find(starMatrix(:, starCol));

while ~isempty(starRow)
    
    % unstar the starred zero
    newStarMatrix(starRow, starCol) = 0;
    
    % find primed zero in row
    primeRow = starRow;
    primeCol = find(primeMatrix(primeRow, :));
    
    % star the primed zero
    newStarMatrix(primeRow, primeCol) = 1;
    
    % find starred zero in column
    starCol = primeCol;
    starRow = find(starMatrix(:, starCol));
    
end
starMatrix = newStarMatrix;

primeMatrix(:) = 0;
coveredRows(:) = 0;

% move to step 2
[assignment, distMatrix, starMatrix, primeMatrix, coveredColumns, coveredRows] = step2__(distMatrix, starMatrix, primeMatrix, coveredColumns, coveredRows, minDim);


% Step 5: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [assignment, distMatrix, starMatrix, primeMatrix, coveredColumns, coveredRows] = step5__(distMatrix, starMatrix, primeMatrix, coveredColumns, coveredRows, minDim)

% find smallest uncovered element
uncoveredRowsIndex    = find(~coveredRows');
uncoveredColumnsIndex = find(~coveredColumns);
[s, index1] = min(distMatrix(uncoveredRowsIndex,uncoveredColumnsIndex));
[s, index2] = min(s); %#ok
h = distMatrix(uncoveredRowsIndex(index1(index2)), uncoveredColumnsIndex(index2));

% add h to each covered row
index = find(coveredRows);
distMatrix(index, :) = distMatrix(index, :) + h;

% subtract h from each uncovered column
distMatrix(:, uncoveredColumnsIndex) = distMatrix(:, uncoveredColumnsIndex) - h;

% move to step 3
[assignment, distMatrix, starMatrix, primeMatrix, coveredColumns, coveredRows] = step3__(distMatrix, starMatrix, primeMatrix, coveredColumns, coveredRows, minDim);


function table = maketbl(y)

maxlevels = max(y);
[counts values] = hist(y,(1:maxlevels));
total = sum(counts);
percents = 100*counts./total;
table = [values' counts' percents'];


function Z = squareform(Y,dir)

if ~(isnumeric(Y) || islogical(Y)) || ndims(Y) > 2
    error('stats:squareform:BadInput','Y must be a numeric or logical vector or matrix.');
end

[m, n] = size(Y);
if nargin<2 || isempty(dir)
    if isvector(Y)
        dir = 'tomatrix';
    else
        dir = 'tovector';
    end
end
okdirs = {'tovector' 'tomatrix'};
j = strmatch(dir,okdirs);
if ~isscalar(j)
    error(message('stats:squareform:BadDirection'));
end
dir = okdirs{j};

switch(dir)
    case 'tomatrix'
        if ~isvector(Y)
            error('stats:squareform:BadInput',...
                'Y must be a vector to convert it to a distance matrix');
        end
        if m~=1
            Y = Y';
            n = m;
        end
        
        m = ceil(sqrt(2*n)); % (1 + sqrt(1+8*n))/2, but works for large n
        if m*(m-1)/2 ~= n
            error('stats:squareform:BadInput',...
                'The size of the vector Y is not correct');
        end
        
        if islogical(Y)
            Z = false(m);
            if m>1
                Z(tril(true(m),-1)) = Y;
                Z = Z | Z';
            end
        else % isnumeric(Y)
            Z = zeros(m,class(Y));
            if m>1
                Z(tril(true(m),-1)) = Y;
                Z = Z + Z';
            end
        end
        
    case 'tovector'
        if m~=n || ~all(diag(Y)==0)
            error('stats:squareform:BadInput',...
                'The distance matrix Z must be square with 0 along the diagonal.');
        end
        
        Z = Y(tril(true(n),-1));
        Z = Z(:)';                 % force to a row vector, even if empty
end

