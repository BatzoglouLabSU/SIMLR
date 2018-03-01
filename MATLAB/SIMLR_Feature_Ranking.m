function [aggR,pval] = SIMLR_Feature_Ranking(A,X)

%%% A is the similarity matrix by SIMLR
%%% X is the data of size nxp

%%%glist is the ranking of genes
%addpath('/Users/Bo/Documents/PHD_in_Stanford/Work_Serafim/SingleCell/GB_similaritylearning/RankAggreg')
for i = 1:100
    i
    index = randperm(length(A));
    index = index(1:round(length(A)*0.9));
    Ai = A(index,index);
    Xi = X(index,:);
    yscore(i,:) = (LaplacianScore(Xi,Ai))';
end

yscore = 1-yscore;
glist = (yscore' - min(yscore(:))+eps)/(max(yscore(:))-min(yscore(:))+eps);
[aggR, pval, rowNames] = aggregateRanks(glist);
[~,aggR] = sort(pval);
pval = pval(aggR);

end


function [aggR, pval, rowNames] = aggregateRanks(R, N, complete, topCutoff)
%AGGREGATERANKS Aggregate ranked lists using traditional and robust methods
%   [aggR, pval, rowNames] = AGGREGATERANKS(R,N,method,complete,topCutoff)
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%   INPUTS:
%       R           -> matrix representation: numeric matrix with as many
%                   rows as the number of unique elements and with as many
%                   columns as the number of ranked lists. All entries have
%                   to be on interval [0,1]. Smaller numbers correspond to
%                   better ranking.
%                   -> list representation: cell array of
%                   cells with strings or vectors with numbers - in this
%                   case R is transformed into numeric rank matrix.
%       N           number of ranked elements, default is the number of
%                   unique elements in R
%       method      rank aggregation method. Could be one of the following:
%                   'min', 'median', 'mean', 'geom.mean', 'stuart', 'RRA'.
%                   Default is 'RRA' (Robust Rank Aggregation).
%       complete    1 - rankings are complete (there is perfect match
%                   between sets of rankings)
%                   0 - default; rankings are incomplete.
%       topCutoff   vector of cutoff values that limit the number of
%                   elements in the input list.
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%   OUTPUTS:
%       aggR        vector of aggregated ranks/scores (it equals pval for
%                   methods 'stuart' and'RRA')
%       pval        p-values (relevant only for 'mean','stuart','RRA')
%       rowNames    if R contains lists, rowNames contains their unique
%                   names in the same order as the values of aggR
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%   EXAMPLES
%       % Lets have three ordered lists of names.
%       R = {   {'Joe', 'Bob', 'Lucy', 'Mary'}, ...
%               {'Bob', 'Joe', 'Lucy', 'Mary'}, ...
%               {'Lucy', 'Mary', 'Bob', 'Joe'}}
%
%       % We can also use numerical vectors instead of cell of strings.
%       % R = { [1,2,3,4], [2,1,3,4], [3,4,2,1] }
%
%       % Obtain aggregated ranks with method 'RRA' (default).
%       [aggR, pval, rowNames] = aggregateRanks(R)
%
%       % Or, equivalently, use explicit parameters definition.
%       [aggR, pval, rowNames] = aggregateRanks(R, [], 'RRA')
%
%       % We can also compute a matrix with ranks first ...
%       rmat = rankMatrix(R)
%       % ... and then pass it to the aggregation method.
%       [aggR, pval, rowNames] = aggregateRanks(rmat)
%
%       % A case of incomplete lists.
%       R = {   {'Joe', 'Bob', 'Lucy', 'Mary'}, ...
%               {'Bob', 'Joe', 'Lucy',       }, ...
%               {'Lucy', 'Mary'              }}
%
%       % Lets compute mean ranks. Mind the fourth parameter, which
%       % indicates completeness of the lists. Note the return values; aggR
%       % contains average across the ranks, while pval contains the
%       % statistical significance (p-value) of mean ranks.
%       [aggR, pval, rowNames] = aggregateRanks(R,[],'mean',0)
%
%       % We can also say that only top k elements are presented in data
%       % by setting the topCutoff to [1,0.75,0.5].
%       [aggR, pval, rowNames] = aggregateRanks(R,[],'RRA',0,[1,0.75,0.5])
%
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%   See also RANKMATRIX.
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%   Copyright (2013) Nejc Ilc <nejc.ilc@gmail.com>
%   Based on R package RobustRankAggreg written by Raivo Kolde.
%   Reference:
%     Kolde, R., Laur, S., Adler, P., & Vilo, J. (2012).
%     Robust rank aggregation for gene list integration and meta-analysis.
%     Bioinformatics, 28(4), 573-580
%   Revision: 1.0 Date: 2013/05/16
%--------------------------------------------------------------------------

rowNames = [];

if ~exist('R','var') || isempty(R)
    error('Input parameter R is missing!');
end

if ~exist('N','var')
    N=[];
end

if ~exist('complete','var') || isempty(complete)
    complete = 0;
end

% Input parameter R determination
if iscell(R)
    [rmat, rowNames] = rankMatrix(R, N, complete);
elseif ismatrix(R)
    if all(max(R,[],1)<=1) && all(min(R,[],1)>0)
        rmat = R;
    else
        error('Columns of matrix R can only contain numbers from interval (0,1].');
    end
else
    error('R should be cell (of lists) or matrix (of ranks).');
end



if ~exist('topCutoff','var') || isempty(topCutoff)
    topCutoff = NaN;
end

pval = NaN;


        aggR = rhoScores(rmat, topCutoff);
        pval = aggR;
        

end


function rho = rhoScores(r, topCutoff)
%RHOSCORES Compute Rho scores for rank vector
%   rho = RHOSCORES(r, topCutoff) 
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -   
%   INPUTS:
%       r           vector of normalized rank values on interval [0,1]
%       topCutoff   a vector of cutoff values used to limit the number of 
%                   elements in the input lists
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%   OUTPUTS:
%       rho         a vector of rho values, corrected against bias from
%                   multiple hypothesis testing (Bonferroni correction).
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%   See also CORRECTBETAPVALUES, THRESHOLDBETASCORE, AGGREGATERANKS.
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%   Copyright (2013) Nejc Ilc <nejc.ilc@gmail.com> 
%   Based on R package RobustRankAggreg written by Raivo Kolde. 
%   Reference:
%     Kolde, R., Laur, S., Adler, P., & Vilo, J. (2012).
%     Robust rank aggregation for gene list integration and meta-analysis.
%     Bioinformatics, 28(4), 573-580
%   
%   Revision: 1.0 Date: 2013/05/16
%--------------------------------------------------------------------------
    if ~exist('topCutoff','var') || isempty(topCutoff)
        topCutoff = NaN;
    end
    
    if isvector(r)
        rows = 1;
        r = r(:)'; % force row vector form
    else
        rows = size(r,1);
    end
    
    rho = nan(rows,1);
    
    for rInd = 1:rows
        r1 = r(rInd,:);
        
        if(isnan(topCutoff(1)))
            x = betaScores(r1);
            % Correct using Bonferroni method.
            rho(rInd) = correctBetaPvalues( min(x), sum(~isnan(x)));
        else            
            r1 = r1(~isnan(r1));
            r1(r1 == 1) = nan;
            % Consider thresholds in topCutoff vector.
            x = thresholdBetaScore(r1,[],[],topCutoff);
            % Correct using Bonferroni method.
            rho(rInd) = correctBetaPvalues( min(x), length(r1));
        end
    end
end

function pval = correctBetaPvalues(p,k)
%CORRECTBETAPVALUES Compute p-values based on Beta distribution
%   pval = CORRECTBETAPVALUES(p,k) 
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%   Copyright (2013) Nejc Ilc <nejc.ilc@gmail.com> 
%   Based on R package RobustRankAggreg written by Raivo Kolde. 
%   Reference:
%     Kolde, R., Laur, S., Adler, P., & Vilo, J. (2012).
%     Robust rank aggregation for gene list integration and meta-analysis.
%     Bioinformatics, 28(4), 573-580
%   
%   Revision: 1.0 Date: 2013/05/16
%--------------------------------------------------------------------------
    pval = betacdf(p,1,k);
end
