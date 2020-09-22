function convertToKalman(filename)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    [filepath,name,ext] = fileparts(filename);
    
    if strcmp(ext,'.czi')
        img = bfLoadTif(filename);
    else strcmp(ext,'.tif')
        img = loadTif(filename);
    end
    
    imgfilt = Kalman_Stack_Filter(single(img));
    [m,n,t] = size(imgfilt);
    imgfilt = reshape(imgfilt,m,n,1,t);
    v = VideoWriter([filepath '\' name '_Kalman.mp4'],'MPEG-4');
    open(v)
    writeVideo(v,mat2gray(imgfilt))
    close(v)

end

