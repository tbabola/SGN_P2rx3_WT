%%%
list = loadFileList('.\Data\*\*.czi')

for i=1:size(list,1)
    [fp,name,ext] = fileparts(list{i});
    imgdata = bfopen([fp '\' name ext]);
    imgdata = imgdata{1,1}; %this is where the imagedata is stored
    img = []; index = 1;
    for i = 1:600%1:size(imgdata,1) 
        img(:,:,index) = imgdata{i,1};
        index = index + 1;
    end
    img = Kalman_Stack_Filter(img);
    img = img ./ max(img,[],'all');
    
    v = VideoWriter([fp '\' name '_Kalman.mp4'],'MPEG-4');
    open(v);
    for i = 1:size(img,3)
        writeVideo(v,img(:,:,i));
    end
    close(v);
    
end
