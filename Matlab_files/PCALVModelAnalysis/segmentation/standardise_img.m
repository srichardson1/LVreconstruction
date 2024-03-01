function imDataN = standardise_img(imData)

imData = double(imData);
Imax = max(imData(:));
Imin = min(imData(:));

imDataN = (imData - Imin)./(Imax - Imin);
imDataN = imDataN.*255;