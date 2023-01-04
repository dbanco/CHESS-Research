function rationals = genRationals(lowLim,hiLim,denLim,numLim)
%genRationals generates a bunuch of rational numbers between two values up to a
% denominator limit
rationals = [lowLim,hiLim];
numVals = size(rationals,2);
while (numVals < numLim)
    N = 2*numVals-1;
    rationals_new = zeros(2,N);
    rationals_new(:,1:2:N) = rationals;
    for i = 2:2:N-1
        if (rationals_new(2,i-1)+rationals_new(2,i+1)) < denLim
            [Nu,De] = rat((rationals_new(1,i-1)+rationals_new(1,i+1))/(rationals_new(2,i-1)+rationals_new(2,i+1)));
            rationals_new(:,i) = [Nu;De];
        end
    end
    % eliminate terms not set
    zeroTerms = rationals_new(1,:) == 0;
    rationals_new(:,zeroTerms) = [];
    
    % reduce fractoins
    fractions = mod(rationals_new(1,:),rationals_new(2,:));
    
    if size(rationals_new,2) == numVals
        break;
    end
    clear rationals
    rationals = rationals_new;
    numVals = size(rationals,2);
    
end



% Remove 1/1 if it's at the beginning still
if rationals(1,1)/rationals(2,1) == 1
    rationals(:,1) = [];
end