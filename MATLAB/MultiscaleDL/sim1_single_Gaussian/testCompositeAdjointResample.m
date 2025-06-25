function testCompositeAdjointResample()
    % Composite op: D3 * U2, Adjoint: U3 * D2
    K = 1;
    scales = cell(K,1);
    scales{1} = genRationals([0;1],[1;1],16,16, 1/8);
    J = size(scales{1},2);
    c = scales{1};

    k = 1;
    count = 0;
    M = 71;
    center = ceil((M+1)/2);
    for N = 117:2:125

        % Create random signals
        % x = randn(1,N,K);
        x = zeros(1,M,1);
        x(1,11:51,1) = [1:20,21,20:-1:1];
        y = randn(1,N,J);

        NormVals = ones(J,1);
        [A_x,NormVals,centers] = reSampleCustomArrayCenter(N,x,scales,center,NormVals);
        A_T_y = reSampleTransCustomArrayCenter(M,y,scales,centers,NormVals);


        % === TEST ===
        lhs = sum(A_x .* y,'all');
        rhs = sum(x .* A_T_y);
    
        if abs(lhs - rhs) > 1e-10
            fprintf(' N = %i \n',N)
            count = count + 1;
        end
        k = k + 1;
    end
    fprintf('%i/%i errors \n',count,k-1)
end
