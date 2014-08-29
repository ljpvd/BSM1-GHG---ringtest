function [Xcf X CXcf errmin]=lhs_imancon(nsample,nvar,xlu,C)
% LHS with correlation method of Iman & Conover
% Iman, R. L., and W. J. Conover. 1982. A Distribution-free Approach to Inducing Rank Correlation
%      Among Input Variables. Communications in Statistics B 11:311-334
%
% Input:
%   C    : correlation matrix of the variables (nvar,nvar)
%   nsample : no. of samples
%   nvar    : no. of variables
%   xlu= lower and upper bound of parameter values
% Output:
%   X       : LHS stratified random samples (nsample,nvar) (units P:0-1)
%   Xcf     : LHS sampling with corretion control (nsample,nvar)(units
%   P:0-1)
%   Gurkan Sin (2011), DTU @ Technical University of Valencia.

% Here we implement the Helton and Davis algorithm of Iman-Conover rank
% transformation based correlation-induction method.

% First step: generate nsample of van der Waarden scores.
% Note that it was found out that it provides better results if LHS is used
% to generate the S matrix instead of van der Waarden scores.
errmin = 1000 ; % this is arbitrarily chosen to iterate until we got a good convergence in the correlation matrix.
% Generate LHS sample with no correlation
XS = lhsdesign(nsample,nvar);
% From probability to values via inverse CDF
lb=xlu(1,:)'*ones(1,nsample);
ub=xlu(2,:)'*ones(1,nsample);
X = unifinv(XS,lb',ub');  % this part should be edited by the user depending on the prior prob distribution.

%then run iteration to impose correlation control by C
for it=1:100
    for j=1:nvar
        stemp=randperm(nsample)./(nsample+1);
        S(:,j)=norminv(stemp,0,1);       % inverse of standard normal distribution.
    end
    
    % % the below is optional as well as above (Stein and Budiman uses below
    % method
    % xm = zeros(nvar,1);
    % xs = eye(nvar);
    % S = lhsnorm(xm,xs, nsample);
    % %% The correlation matrix should be given as rank-based.
    
    %C = corr2;
    %% Decompose the correlation matrix: C: Note that this must be positive
    %% definite
    % using modified cholesky for corr matrix that is not quite positive definite
    P = chol(C); % Upper triangular matrix
    %%below is an alternative cholesky decomposition
    %     [L,D,E1]=mchol(C);
    %     P=L*sqrt(D);

    %% The correlation matrix of S should be an identity matrix.
    %% This is ensured as follows:
    E = corrcoef(S);
    Q = chol(E); % Upper triangelur matrix
    Q=Q'; % Transposed to lower triangular matrix
  %an alternative cholesky decomposition:
%     [L,D,E1]=mchol(E);
%     Q=L*sqrt(D);
    % finally calculate the transformed matrix
    Str = S*inv(Q)'*P ;
    
    % First calculate the rank of Str's elements
    [srtS ix1] =sort(Str);
    
    rStr = zeros(nsample,nvar);
    for i=1:nvar
        rStr(ix1(:,i),i)=1:nsample;
    end
    
    % in fact the difference for the pearson correlation is zero:
    % CrStr = corrcoef(rStr);
    % CS = corrcoef(Str);
    % err1 = CS - C ;
    
    %% Now we need to rearrange the elements of X such that it corresponds with
    %% the rank order of the elements of Str, i.e. rStr:
    
    [srtX ix2] =sort(X);
    
    for i=1:nvar
        Xcorr(:,i)=srtX(rStr(:,i),i);
    end
    
    CX = corrcoef(Xcorr);
    err2 = sum(abs(CX - C)) ;
    
    % update
    if(err2<errmin)
        Xcf=Xcorr;
        CXcf = CX;
        errmin = err2;
    end
    
end
