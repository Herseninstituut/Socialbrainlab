%
function filter3image

Miji
IJ = ij.IJ;

ori=double(MIJ.getCurrentImage());
for i=1:size(ori,3)-2
    flt(:,:,i)=round((ori(:,:,i)+ori(:,:,i+1)+ori(:,:,i+2))/3);
end
MIJ.createImage(single(flt));