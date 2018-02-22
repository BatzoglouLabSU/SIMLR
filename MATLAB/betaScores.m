function p = betaScores(r)
%BETASCORES Compute Beta scores for rank vector
%   p = BETASCORES(r) 
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -   
%   INPUTS:
%       r   vector of normalized rank values on interval [0,1]
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%   OUTPUTS:
%       p   a vector of p-values that corresponds to the sorted input 
%           vector. The NaN-s are moved to the end.
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
    n = sum(~isnan(r));
    p = nan(1,length(r));
    % Sort the values.
    r = sort(r);
    % Get the order statistics and calculates p-values for each of the
    % order statistics. These are based on their expected distribution
    % under the null hypothesis of uniform distribution.
    p(1:n) = betacdf(r(1:n),1:n,n:-1:1);
end
