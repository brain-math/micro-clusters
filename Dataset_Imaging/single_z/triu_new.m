function x = triu_new(x, needDiag, NonNaNElementsOnly)

if nargin == 1
    error('input x, needDiag, and NonNaNElementsOnly.');
end

% Constants.
m = cast(size(x,1),coder.internal.indexIntClass);
n = cast(size(x,2),coder.internal.indexIntClass);
ONE = ones(coder.internal.indexIntClass);

k = 0;

% Calculate limits.
if isempty(x) || 1 - k >= m
    % Trivial case.  No change to x.
    return
elseif k > 0
    istart = ONE;
    if k < n % Guarantee 1 <= jstart <= n.
        jstart = cast(k,coder.internal.indexIntClass);
    else
        jstart = n;
    end
else
    istart = cast(2-k,coder.internal.indexIntClass); %  2 <= istart <= m is guaranteed.
    jstart = ONE;
end
nrowsm1 = coder.internal.indexMinus(m,istart);
ncolsm1 = coder.internal.indexMinus(n,jstart);
if nrowsm1 < ncolsm1
    jend = coder.internal.indexPlus(jstart,nrowsm1);
else
    jend = coder.internal.indexPlus(jstart,ncolsm1);
end

% Zero matrix elements.
if needDiag == 0
for j = ONE:jend+1
    for i = istart-1:m
        x(i,j) = NaN;
    end
    if j >= jstart
        istart = coder.internal.indexPlus(istart,ONE);
    end
end
elseif needDiag == 1
for j = ONE:jend
    for i = istart:m
        x(i,j) = NaN;
    end
    if j >= jstart
        istart = coder.internal.indexPlus(istart,ONE);
    end
end
else
error('needDiag should be 0 or 1');
end

if NonNaNElementsOnly == 1
    x = x(:);
    x = x(~isnan(x));
end
