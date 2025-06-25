function testCompositeAdjoint2()
    % Composite op: D3 * U2, Adjoint: U3 * D2
    K = 1;
    scales = cell(K,1);
    scales{1} = genRationals([0;1],[1;1],16,16, 1/8);
    J = size(scales{1},2);
    c = scales{1};

    k = 1;
    count = 0;
    M = 71;
    for N = 117:2:1256
        for i = 1:size(c,2)
            c1 = c(1,i); c2 = c(2,i);
            center = ceil((N+1)/2);
        
            % Create random signals
            d = zeros(1,M,1);
            d(1,11:51,1) = [1:20,21,20:-1:1];
            y = randn(1,N);

            % Forward op: D3 * U2
            [u2, center_up] = uSampleCenter2(d, c1, c1*center);  % len_u2
            [A_x, center_dwn] = dSampleCenter2(u2, c2, center_up);   % N
            Nout = numel(A_x);
            A_x = padarray(A_x,[0,N-Nout],0,'post');
            % y = randn(1, Nout);
            % === ADJOINT PATH ===
            % U3 * D2: y → d2 → x_hat
        
            [d2, center_up2] = uSampleCenter2(y, c2, c2*center_dwn);  % len_u2
            [A_T_y, center_dwn2] = dSampleCenter2(d2, c1, center_up2); % len_x
            % Ndwn2 = numel(A_T_y);
            % Nout2 = round(N);
            shift = center_dwn2-center;
            A_T_y = A_T_y(:,1+shift:M+shift);
        
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
end
