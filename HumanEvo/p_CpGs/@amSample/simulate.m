function simu_as = simulate (obj,mms)

% SIMULATE creates simulations of methylation values
% -------------------------------------------------------------------
% Description:  simulate C and T values based on the parameters of the
%               sample and the methylation map of the reference, and
%               converts them into methylation values.
% Inputs:           {obj} an amSample object with counts data
%                   {mms} an object that contains the reference 
%                   methylation map for the simulations.
%                   
%
% Output:       {TA_convPerc} cell array over chromosomes, containing the
%                   methylation value of each position along each
%                   chromosome.

% © Liran Carmel & Yoav Mathov
% Classification: RoAM
% Last revision date: 29-Jan-2019

simu_as = amSample;

MethMapRef = mms.getmethylation;

degrad_rate = obj.drate.rate.global;
simu_as = simu_as.set ('No_Cs',obj.No_Cs);

for chr = 1:length (obj.No_Cs)
    totCT = obj.No_Cs{chr} + obj.No_Ts{chr};
    idx = find(totCT > 0);
    temp = binornd(double(totCT(idx)),...
        0.01*degrad_rate*MethMapRef{chr}(idx));
    simu_as.No_Ts{chr} = zeros(length(MethMapRef{chr}),1);
    simu_as.No_Ts{chr}(idx) = temp;
end
