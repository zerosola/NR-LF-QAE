function [EPI_Gabor_width, EPI_Gabor_height] = EPI_Gabor(pos, lf_num, rows, cols)
for k = 1:lf_num
    k
    for i = 1:9
        for j = 1:9
            im = imread(strcat(pos,'\',num2str(k),'\',num2str(i),num2str(j),'.bmp'));%SAI
            dis_lf(:,:,j) = double(rgb2gray(im));
        end
        Slice_epi_width = permute(dis_lf,[1,3,2]);
        % imshow(squeeze(Slice_epi_width(200,:,:)),[])
        for numi = 1 : rows
            [col_mean,col_entropy,col_skewness,col_kurtosis] = lpc_si(squeeze(Slice_epi_width(numi,:,:)));
            EPI_Gabor_width(k,numi,1:4*3) = [col_mean,col_entropy,col_skewness,col_kurtosis];
        end
    end
    
    dis_lf = [];
    for j = 1:9
        for i = 1:9
            im = imread(strcat(pos,'\', num2str(k),'\',num2str(i),num2str(j),'.bmp'));%SAI
            dis_lf(:,:,i) = double(rgb2gray(im));
        end
        Slice_epi_height = permute(dis_lf,[2,3,1]);
%         imshow(squeeze(Slice_epi_height(200,:,:)),[])
        for numi = 1 : cols
            [col_mean,col_entropy,col_skewness,col_kurtosis] = lpc_si(squeeze(Slice_epi_height(numi,:,:)));
            EPI_Gabor_height(k,numi,1:4*3) = [col_mean,col_entropy,col_skewness,col_kurtosis];
        end
    end
    
end

save EPI_Gabor_width.mat EPI_Gabor_width
save EPI_Gabor_height.mat EPI_Gabor_height