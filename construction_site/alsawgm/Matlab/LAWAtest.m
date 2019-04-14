addpath(genpath(getenv('MATLAB_APPROXIMATIONTB')))
clear all

d       = 4;
r       = 3;
rt      = 2;

Alawa   = cell(d, 1);
for i = 1:d
    Alawa{i} = ones(r, r);
end

blawa   = cell(d, 1);
for i = 1:d
    blawa{i} = ones(r, 1);
end

x0lawa    = cell(d-1, 1);
x0lawa{1} = ones(r*rt, r);
x0lawa{2} = ones(r*rt, rt);
x0lawa{3} = ones(r, rt);

r0lawa       = [rt, rt, 1];

maxIt        = 5;
tol          = 1e-02;
stag         = 1e-02;

% Test conversions
A  = LAWATENSALGConvertor.lawa2tensalgLaplaceOperator(Alawa);
b  = LAWATENSALGConvertor.lawa2tensalgCanonicalTensor(blawa);
x0 = LAWATENSALGConvertor.lawa2tensalgHTCores2TT(x0lawa, r0lawa);

% Apply solver
s       = TTTensorALSLinearSolver(OperatorTTTensor(A),...
                                  TTTensor(b),...
                                  'x0', x0,...
                                  'maxIterations', maxIt,...
                                  'tolerance', tol,...
                                  'stagnation', stag,...
                                  'algorithm', 'dmrg');
x       = s.solve();
