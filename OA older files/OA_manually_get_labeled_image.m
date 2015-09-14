function [Ifin2]=OA_manually_get_labeled_image(Itmp2,I,IG,IY) 

%now also displays the GFP

[x_size,y_size]=size(Itmp2);
Itmp2new=Itmp2;
 Ifinal=zeros(x_size,y_size);
 checkbutton=1;
 
 figure(2);imshow(I);colorbar
 figure(5);imshow(IY);colorbar
 figure(6);imshow(IG);colorbar
 while checkbutton~=9
     figure(1);imagesc(Itmp2new);colorbar
     imtext=['left button -> remove, right button add, wheel for approve/disapprove, click upper right when finished'];
     text(10,10,imtext,'color',[0 1 1])
     [x1,y1,button] = ginput(1);
     if x1>(0.9.*y_size) && y1<(0.1.*x_size) %upper right corner
     checkbutton=9;
     else
         checkbutton=button;
     end
     Itmphere=zeros(x_size,y_size);    
     Itmphere(round(y1),round(x1))=1;
         intr_area=bwmorph(Itmphere,'dilate',1.5);%was2
     if button ==1 %cut two cells
     Itmp2new=(  (Itmp2new.*double(~intr_area))+0.*double((double(intr_area).*double(Itmp2new+1))>0));
     elseif button ==2 %remove the marked cell and find its "correct dimension"
         tmpL2=bwlabel(Ifinal);
         if tmpL2(round(y1),round(x1))>0 %we marked a moved obj
             onum=tmpL2(round(y1),round(x1));
             Itmp2new=Itmp2new+double(tmpL2==onum).*Itmp2;
             Ifinal=Ifinal-(tmpL2==onum);
         else %take the obj away
         tmpL=bwlabel(Itmp2new>0);
         b=unique(tmpL.*Itmphere);
         Ifinal=Ifinal+(tmpL==b(2));
         figure(3);imshow(Ifinal);
         Itmp2new=( Itmp2new.*double(~(tmpL==b(2))) );
         figure(2);imshow(I.*uint8(~Ifinal));colorbar
         end
     elseif button==3 %fill in line      
     Itmp2new=(  (Itmp2new.*double(~intr_area))+6.*double((double(intr_area).*double(Itmp2new+1))>0));
     end
 end
 
 Lf=bwlabel(Ifinal,4);
 
 Ifin2=zeros(x_size,y_size);    
 for i=1:max(max(Lf))
     Ifin2=Ifin2+imfill(bwmorph(Lf==i,'close',1),'holes');
 end
 %figure(4);imshow(Ifin2)