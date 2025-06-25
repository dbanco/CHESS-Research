
% Composite op: D3 * U2, Adjoint: U3 * D2
K = 1;
scales = cell(K,1);
scales{1} = genRationals([0;1],[1;1],16,16, 1/8);
J = size(scales{1},2);
c = scales{1};

figure
hold on
k = 1;
count = 0;
M = 71;
for N = 117
    for i = 1:size(c,2)
        scales_i = {c(:,i)};
        c1 = c(1,i);
        c2 = c(2,i);
    
        % Create random signals
        d = zeros(1,M,1);
        center = ceil((M+1)/2);
        d(1,16:56,1) = [1:20,21,20:-1:1];
        y = randn(1,N);
        

        [A_x,NormVals,Shifts] = reSampleCustomArrayCenter3(N,d,scales_i,center); 
        A_T_y = reSampleTransCustomArrayCenter3(M,y,scales_i,center,NormVals,Shifts);
        plot(A_x)
        % === TEST ===
        lhs = sum(A_x .* y);
        rhs = sum(d .* A_T_y);
    
        % fprintf('Composite adjoint test:\n');
        % fprintf('lhs = %.12e, rhs = %.12e, error = %.2e\n', ...
        %     lhs, rhs, abs(lhs - rhs));
        
        if abs(lhs - rhs) > 1e-10
            fprintf('c1 = %i, c2 = %i, N = %i \n',c1,c2,N)
            count = count + 1;
        end
        k = k + 1;
    end
end
fprintf('%i/%i errors \n',count,k-1)

