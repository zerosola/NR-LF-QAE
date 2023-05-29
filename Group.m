function [svd_proportion] = Group(pos, lf_num)

mesh_idx(1,:) = [41, 32, 23, 14, 5];
mesh_top_left = mesh_idx(1,:);
mesh_idx_temp = mesh_idx(1,:);
for j = 1 :4
    mesh_idx_temp = flip(mesh_idx_temp);
    mesh_top_left =[mesh_top_left, mesh_idx_temp-j];
end

mesh_top_right = mesh_idx(1,:);
mesh_idx_temp = mesh_idx(1,:);
for j = 1 :4
    mesh_idx_temp = flip(mesh_idx_temp);
    mesh_top_right =[mesh_top_right, mesh_idx_temp+j];
end

mesh_idx(2,:) = [41, 50, 59, 68, 77];
mesh_down_left = mesh_idx(2,:);
mesh_idx_temp = mesh_idx(2,:);
for j = 1 :4
    mesh_idx_temp = flip(mesh_idx_temp);
    mesh_down_left =[mesh_down_left, mesh_idx_temp-j];
end

mesh_down_right = mesh_idx(2,:);
mesh_idx_temp = mesh_idx(2,:);
for j = 1 :4
    mesh_idx_temp = flip(mesh_idx_temp);
    mesh_down_right =[mesh_down_right, mesh_idx_temp+j];
end

mesh_num = [mesh_top_left;  mesh_top_right;  mesh_down_left;  mesh_down_right];



for k = 1:lf_num
    k
    for i = 1:9
        for j = 1:9
            im = imread(strcat(pos,'\',num2str(k),'\',num2str(i),num2str(j),'.bmp'));%SAI
            ALL_LF(:,:,(i-1)*9+j) = double(rgb2gray(im));
        end
    end
    
    Block_Size = 8;
    Block_Size_Stride = 8;
    Search_Size = 7;
    Search_Stride = (Search_Size-1)/2;

    N     =  size(ALL_LF,1);
    M     =  size(ALL_LF,2);
    D     =  size(ALL_LF,3);
    
    
    for zi = [Block_Size : Block_Size_Stride : N]
        for zj = [Block_Size : Block_Size_Stride : M]
            Block_Ref_LF = ALL_LF(zi-Block_Size+1:zi,zj-Block_Size+1:zj, 41);
            error_distance = zeros(Search_Size,Search_Size);
            
            for meshn = 1 : 4
                delx = 0; dely = 0;
                for j = mesh_num(meshn,:)
                    
                    for Searchi = -Search_Stride:Search_Stride
                        for Searchj = -Search_Stride:Search_Stride
                            search_xmax = min(max(zi+Searchi+delx,Block_Size),N);
                            search_ymax = min(max(zj+Searchj+dely,Block_Size),M);
                            Search_Block = ALL_LF(search_xmax-Block_Size+1:search_xmax,search_ymax-Block_Size+1:search_ymax,  j);
                            error_distance(Searchi+(Search_Size+1)/2,Searchj+(Search_Size+1)/2) = mean(abs(Block_Ref_LF(:)-Search_Block(:)));
                        end
                    end
                    [val,ind] = sort(error_distance(:));
                    [x,y]=ind2sub([Search_Size,Search_Size],ind(1));
                    search_xmax = min(max(zi+x- (Search_Size+1)/2 +delx,Block_Size),N);
                    search_ymax = min(max(zj+y- (Search_Size+1)/2 +dely,Block_Size),M);
                    Search_Block = ALL_LF(search_xmax-Block_Size+1:search_xmax,  search_ymax-Block_Size+1:search_ymax,  j);
                    Match_Patch(zi-Block_Size+1:zi,zj-Block_Size+1:zj,j)= Search_Block;
                    delx = delx + x - (Search_Size+1)/2 ;
                    dely = dely + y - (Search_Size+1)/2 ;
                    
                end 
            end
            
        end
    end
    
    svd_count = 0;
    Block_Size_Area = Block_Size*Block_Size;
    for zi = [Block_Size : Block_Size_Stride : N]
        for zj = [Block_Size : Block_Size_Stride : M]
            [~, SG_V, ~] = svd(reshape(Match_Patch(zi-Block_Size+1:zi,zj-Block_Size+1:zj,:), Block_Size_Area, [] ));            
            SG_V = diag(SG_V);
            svd_count = svd_count + 1;
            svd_proportion(k, svd_count) = SG_V(1)/sum(SG_V);
        end
    end

end

save svd_proportion.mat svd_proportion

