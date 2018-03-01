clear; clc
G1 = im2double(imread('images/part2/apple.jpg'));
figure; imshow(G1)
laplace_filter = fspecial('laplacian', 0.5);
laplace_filter = [0, -1, 0; -1, 4, -1; 0, -1, 0];
G1_laplace = imfilter(G1, laplace_filter, 'replicate'); 
G1_laplace = G1_laplace + 0.5;
figure; imshow(G1_laplace)




