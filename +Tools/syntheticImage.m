function img = syntheticImage(m, n)
%SYNTHETIC_IMAGE Create a deterministic grayscale test image.
%   IMG = SYNTHETIC_IMAGE(M, N) returns an M-by-N image in [0,1]. It is used
%   by create_tv_problem when MATLAB's cameraman.tif is not available.
%
%   Input:
%       m,n - image height and width.
%
%   Output:
%       img - M-by-N double array with simple geometric features.
[X, Y] = meshgrid(linspace(-1, 1, n), linspace(-1, 1, m));
img = 0.18 + 0.22 * (X > -0.7 & X < 0.65 & Y > -0.85 & Y < 0.75);
img = img + 0.35 * ((X + 0.25).^2 + (Y + 0.1).^2 < 0.22);
img = img + 0.25 * (abs(X - 0.45) < 0.08 & Y > -0.65 & Y < 0.45);
img = img + 0.16 * (sin(18*X) > 0.75 & Y < -0.35);
img = min(max(img, 0), 1);
end
