No_As = cell(1,22);
No_Cs = cell(1,22);
No_Gs = cell(1,22);
No_Ts = cell(1,22);
for chr = 1:22
    vec = as.No_As{chr};
    vec = vec(1:2:end);
    No_As{chr} = vec;
    vec = as.No_Cs{chr};
    vec = vec(1:2:end);
    No_Cs{chr} = vec;
    vec = as.No_Gs{chr};
    vec = vec(1:2:end);
    No_Gs{chr} = vec;
    vec = as.No_Ts{chr};
    vec = vec(1:2:end);
    No_Ts{chr} = vec;
end
as = set(as,'No_As',No_As,'No_Cs',No_Cs,'No_Gs',No_Gs,'No_Ts',No_Ts,...
    'coord_per_position',ones(1,22));