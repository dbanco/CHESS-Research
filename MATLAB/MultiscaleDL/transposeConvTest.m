

K = 1;
U = 3;
x = [1 2;3 4];
A = [1 2;3 4];
Apad = [0 0 0;0 1 2;0 3 4];
xx = [1;2;3;4];
AA = [0 0 0 1;0 0 1 2;0 1 0 3; 1 2 3 4];

% Ax and A*x
y1 = convn(Apad,rot90(x,2),'valid');
y2 = AA*xx;

% Atx and transpose of A*x
yt1 = convn(padarray(convn(Apad,rot90(x,2),'valid'),[1 1 1],0,'pre'),rot90(x,2),'valid');
yt2 = AA'*AA*xx;





Ax = zeros(3,3,3);
Ax(1,1,1) = 1;
Ax(2,2,2) = 1;
Ax(3,3,3) = 1;
u = ones(3,3,3);
v = ones(3,3,3);
AxPad = padarray(Ax,[1 1 1],0,'pre');
[Fx,Fy,Ft] = diffDict(Ax,K,U);
[uFx,vFy,Ftt] = diffDictTrans(Ax,K,U,u,v);
out1 = opticalFlowOp(Ax,u,v,K);
out2 = opticalFlowOp(Ax,u,v,K,1);

norm(Fx(:)-uFx(:))
norm(Fy(:)-vFy(:))
norm(Ftt(:)-Ft(:))

%% Checking transposes again
X = zeros(3,3);
X(1,1) = 1;
X(2,2) = 0;
Psi = [  1  0 0  0  0 0  0  0 0;
        -1  1 0  0  0 0  0  0 0;
         0 -1 1  0  0 0  0  0 0;
         1  0 0  1  0 0  0  0 0;
        -1  1 0 -1  1 0  0  0 0;
         0 -1 1  0 -1 1  0  0 0;
         0  0 0  1  0 0  1  0 0;
         0  0 0 -1  1 0 -1  1 0;
         0  0 0  0 -1 1  0 -1 1];

reshape(-Psi*X(:),[3 3])
reshape(-Psi'*X(:),[3 3])


dataPad = padarray(X,[1 1 1],0,'pre');
diffyHS(dataPad)
dataPad = padarray(-X,[1 1 1],0,'post');
diffyHS(dataPad)




