

function[im] = readim(imin)
t = Tiff(imin,'r');
for i3 = 1:8
    t.setDirectory(i3)
    im(:,:,i3) = t.read();
end
close
end
