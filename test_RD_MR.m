function [SNR_out] = test_RD_MR(...
    modes, IFs, clwin, sigma_s, Nfft, smooth_p, SNR_IN, NRep)

N_SNR = length(SNR_IN);
N_sp = length(smooth_p);

[Nr, L] = size(modes);
s_in = sum(modes, 1);

Sub = struct('Cl', zeros(Nr, N_SNR),...
    'MB', zeros(Nr, N_SNR),...
    'New', zeros(Nr, N_SNR, N_sp));
SNR_out.MR = Sub;
SNR_out.LCR = Sub;
SNR_out.RD = Sub;

for n=1:length(SNR_IN)
    for k=1:NRep
        fprintf('snr %d/%d, rep %d/%d\n', n, length(SNR_IN), k, NRep);
        
        noise = randn(L, 1)+1i*randn(L, 1);
        s_noise = sigmerge(transpose(s_in), noise, SNR_IN(n));
        [g, Lh] = create_gaussian_window(L, Nfft, sigma_s);
        X_win = 2*Lh:(L-2*Lh);
        [STFT, omega, ~, QM, ~, tau] = FM_operators(s_noise, L, Nfft, g, Lh, sigma_s);

%         fprintf('Classic, ');
        [Cs_simple] = exridge_mult(STFT, Nr, 0, 0, clwin);
%         Spl_Cl = struct('spline', cell(1, Nr));
%         for m=1:Nr
%             Spl_Cl(m).spline = spline((0:L-1)/L, (Cs_simple(m, :) - 1)*L/Nfft);
%         end
%         [m_SR_Cl, m_LCR_Cl, IF_Cl] = R1_MR_and_LCR_spl(STFT, Spl_Cl, g, Lh, sigma_s, Nr, Nfft, L);
        [m_SR_Cl, m_LCR_Cl, IF_Cl] = R1_MR_and_LCR_grid(STFT, QM, Cs_simple, g, Lh, sigma_s, Nr, Nfft, L);

%         fprintf('VFB MB, ');
        [Cs_VFB_MB] = VFB_MB_exridge_MCS(STFT, sigma_s, QM, 2, Nr);
%         Spl_MB = struct('spline', cell(1, Nr));
%         for m=1:Nr
%             Spl_MB(m).spline = spline((0:L-1)/L, (Cs_VFB_MB(m, :) - 1)*L/Nfft);
%         end
%         [m_SR_MB, m_LCR_MB, IF_MB] = R1_MR_and_LCR_spl(STFT, Spl_MB, g, Lh, sigma_s, Nr, Nfft, L);
        [m_SR_MB, m_LCR_MB, IF_MB] = R1_MR_and_LCR_grid(STFT, QM, Cs_VFB_MB, g, Lh, sigma_s, Nr, Nfft, L);
        
        %% Classic and MB SNR
        for m = 1:Nr
            ref_mode = modes(m, X_win);
            ref_IF = IFs(m, X_win);
            
            x_MR_Cl = snr(ref_mode, m_SR_Cl(m, X_win) - ref_mode);
            x_MR_MB = snr(ref_mode, m_SR_MB(m, X_win) - ref_mode);
            x_LCR_Cl = snr(ref_mode, m_LCR_Cl(m, X_win) - ref_mode);
            x_LCR_MB = snr(ref_mode, m_LCR_MB(m, X_win) - ref_mode);
            x_Cl_RD = snr(ref_IF, IF_Cl(m, X_win) - ref_IF);
            x_MB_RD = snr(ref_IF, IF_MB(m, X_win) - ref_IF);
        
            SNR_out.MR.Cl(m, n) = SNR_out.MR.Cl(m, n) + x_MR_Cl/NRep;
            SNR_out.MR.MB(m, n) = SNR_out.MR.MB(m, n) + x_MR_MB/NRep;
            SNR_out.LCR.Cl(m, n) = SNR_out.LCR.Cl(m, n) + x_LCR_Cl/NRep;
            SNR_out.LCR.MB(m, n) = SNR_out.LCR.MB(m, n) + x_LCR_MB/NRep;
            SNR_out.RD.Cl(m, n) = SNR_out.RD.Cl(m, n) + x_Cl_RD/NRep;
            SNR_out.RD.MB(m, n) = SNR_out.RD.MB(m, n) + x_MB_RD/NRep;
        end

        %% RRP
%         fprintf('RRP RD\n');
        for ns=1:N_sp
            [Spl, ~] = R1_RRP_RD(STFT, QM, omega, tau, L, Nfft, Nr, sigma_s, smooth_p(ns));
            [m_SR_New, m_LCR_New, IF_New] = R1_MR_and_LCR_spl(STFT, Spl, g, Lh, sigma_s, Nr, Nfft, L);
            
            for m = 1:Nr
                ref_mode = modes(m, X_win);
                ref_IF = IFs(m, X_win);
                
                x_MR_New = snr(ref_mode, m_SR_New(m, X_win) - ref_mode);
                x_LCR_New = snr(ref_mode, m_LCR_New(m, X_win) - ref_mode);
                x_RD_New = snr(ref_IF, IF_New(m, X_win) - ref_IF);
                if isnan(x_LCR_New)
                    fprintf("NaN\n");
                end
                SNR_out.MR.New(m, n, ns) = SNR_out.MR.New(m, n, ns) + x_MR_New/NRep;
                SNR_out.LCR.New(m, n, ns) = SNR_out.LCR.New(m, n, ns) + x_LCR_New/NRep;
                SNR_out.RD.New(m, n, ns) = SNR_out.RD.New(m, n, ns) + x_RD_New/NRep;
            end
        end
    end
end

end

