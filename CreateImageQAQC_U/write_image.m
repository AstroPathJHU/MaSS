function write_image(iname,im,Image)
T = Tiff(iname,'w');
T.setTag(Image.ds);
write(T,im);
writeDirectory(T)
close(T)
end
