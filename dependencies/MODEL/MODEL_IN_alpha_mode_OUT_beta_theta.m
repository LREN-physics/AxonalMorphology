function F = MODEL_IN_alpha_mode_OUT_beta_theta(x,alpha,M,G,V,p)
% Mathematical formulation of our proposed model to estimate microscopic
% morphologic features of the axons composing white matter tracts.
% Gets called in EstimateMicrostructure_GroupDelay.m under a non-linear
% least square minimization problem
%
%   INPUTS:
%       alpha (number)      - parameter alpha in g(r)=beta*radius^alpha
%       mode (number)       - mode/peak value of the axon radius distribution (micrometers)
%
%   OUPUTS:
%       beta (number)       - parameter beta in g(r)=beta*radius^alpha
%       theta (number)      - parameter that captures the right-tail of the
%           of the axon radius distribution (micrometers)
%  
% Author: Rita Oliveira
% Email: rita.oliveira.uni@gmail.com
% Laboratory for Research in Neuroimaging (LREN)
% Department of Clinical Neuroscience, Lausanne  University Hospital and University of Lausanne
% Mont-Paisible 16, CH-1011 Lausanne, Switzerland
%
% Last updated: 17/08/2021
%
%------------- BEGIN CODE --------------

% Variables
beta    = x(1);
theta   = x(2);

F       = [];

% Velocity equations (add as much as the input has)
for i=1:length(V)
    v = V(i);
    F = [F;
        -v + ((2*p)/beta)*theta^(1-alpha) ...
        * gamma(M/theta+1-alpha)/gamma(M/theta+1) ...
        * (M/theta+1-alpha)];
end

% G ratio equations (add as much as the input has)
for i=1:length(G)
    g = G(i);
    F  = [F;
        -g^2 + beta^2*theta^(2*alpha) ...
        * gamma(M/theta+1)/gamma(M/theta+1-2*alpha) ... 
        * (M + 3*theta + 2*(theta^2)/M) ...
        / (M + (3 - 4*alpha)*theta + (2+4*alpha^2-6*alpha)*(theta^2)/M)];
end

end
