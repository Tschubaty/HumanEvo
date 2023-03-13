function [no_ti, no_cti] = smooth(obj,chr,W)

% SMOOTH provides smoothed CT and T vectors
% -----------------------------------------
% [no_ti, no_cti] = obj.smooth(chr,winsize)
% -----------------------------------------
% Description:  provides smoothed CT and T vectors using a sliding window.
% Input:        {chr} index of chromosome.
%               {winsize} window size.
% Output:       {no_ti} smoothed T's vector from the chromosome.
%               {no_cti} smoothed CT's vector from the chromosome.

% (c) Liran Carmel
% Classification: SET/GET functions
% Last revision date: 08-Sep-2018

% get T's and CT's vectors
no_ti = obj.getNo_Ts(chr);
no_cti = no_ti + obj.getNo_Cs(chr);

% make the smoothing
tmplt = ones(1,W);
no_ti = nanconv(no_ti,tmplt,'same');
no_cti = nanconv(no_cti,tmplt,'same');