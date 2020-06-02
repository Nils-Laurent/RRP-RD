function [] = compute_GHE()

% RRP_in = RRP_merged;
% RRP_E_in = RRP_E_merged;
% N_RRP = Nm;
% 
% fprintf("compute RRP with global highest energy\n");
% RRP_HE = 1:N_RRP;
% for n=1:L
%     [~, r_sort_vec] = sort(RRP_E_in , 'descend');
%     p = 1;
%     for r=r_sort_vec(Nr+1:end)
%         if RRP_in(r, n) > 0
%             if p > Nr
%                 RRP_HE(r) = 0;
%             end
%             
%             p = p + 1;
%         end
%     end
% end

R_HE = zeros(Nfft, L);
% for r=RRP_HE
%     if r == 0
%         continue;
%     end
%     
%     kr_vec = RRP_in(r, :);
%     for n=1:L
%         if kr_vec(n) == 0
%             continue;
%         end
%         R_HE(kr_vec(n), n) = RRP_E_in(r);
%     end
% end
% 
% figure;
% imagesc(1:L, 1:Nfft, R_HE);
% set(gca,'ydir','normal');
% axis square
% colormap(flipud(gray));
% title("global highest energy");

end

