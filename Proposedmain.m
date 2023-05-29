clear
clc
pos = 'Win5-LID\Distorted\Real';
lf_num = 1;
rows = 434;
cols = 625; 
volume = 81;

tic
parfor z = 1 : 4
    switch(z)
        case 1
            [EPI_Gabor_width, EPI_Gabor_height] = EPI_Gabor(pos, lf_num, rows, cols);
        case 2
            [svd_proportion] = Group(pos, lf_num);
        case 3
            [LBP_features_width, LBP_features_height] = LBP_TOP(pos, lf_num);
        case 4
            [ThreeGabor_features1,ThreeGabor_features3] = ThreeLogGabor(pos, lf_num, rows, cols, volume);
    end
end
toc



