function Jreg = compute_filter_regularizer_penalty(kernel,X,K,J)

Jreg = 0;
for k = 1:K
    i1 = 1 + (k-1)*J;
    i2 = J + (k-1)*J;
    Jreg = Jreg + norm(convn(squeeze(X(1,:,i1:i2,:)),kernel,'same'),'fro');        
end

end