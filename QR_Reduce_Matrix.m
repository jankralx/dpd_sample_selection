function [Xsub, idx] = QR_Reduce_Matrix(X, N)
% [Xsub, idx] = QR_REDUCE_MATRIX(X, N) returns a matrix Xsub with N most 
% significant columns selected by QR decomposition
%
% This function is re-implementation of original Tomas Gotthans method.
%
%   Inputs:
%   =======
%   X  - matrix to be reduced
%
%   N  - number of columns of X to be selected
%
%   Returns:
%   ========
%
%   Xsub  - matrix with N most significant columns
%
%   idx - indexes of selected columns of the original matrix X

% Authors: Tomas Gotthans, Jan Kral <kral.j@lit.cz>
% Date: 6.6.2019

if N >= size(X,2)
    Xsub = X;
    idx = 1:size(X,2);
    return;
end

if ~nnz(X) % X has only zeros and hence no independent columns
    Xsub=[];
    idx=[];
    return
end

% perform QR decomposition
[~, ~, E] = qr(X,0); 

% select N most significant columns
idx = sort(E(1:N));
Xsub = X(:,idx);

end

% original function from Tomas Gotthans

% function [Xsub,idx, diagr]=QR_Reduce_Matrix(X,tol)
% 
%      if ~nnz(X) %X has no non-zeros and hence no independent columns
%          Xsub=[]; idx=[];
%          return
%      end
%      if nargin<2, tol=1e-10; end
%        [Q, R, E] = qr(X,0); 
%        if ~isvector(R)
%         diagr = abs(diag(R));
%        else
%         diagr = R(1);   
%        end
%        %Rank estimation
%        r = find(diagr >= tol*diagr(1), 1, 'last'); %rank estimation
%        idx=sort(E(1:r));
%        Xsub=X(:,idx);