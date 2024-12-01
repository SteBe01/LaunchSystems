function [D, T, T_tot] = maxDiameter(N, I_sp)
    data = importdata("maxD.mat");
    data_filtered = data([data.num_engines]==N);
    D = interp1([data_filtered.Isp], [data_filtered.D_ext], I_sp);
    T = interp1([data_filtered.D_ext], [data_filtered.T_single_engine], D);
    T_tot = T*N;
end