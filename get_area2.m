
%function [vol_c,area_c,maj_appr_axis,min_appr_axis]=get_area2(all_obj,numbM,cell_number)
function [vol_c,area_c,maj_appr_axis,min_appr_axis]=get_area2(new_c_Image2,t_here,cell_number)

%program for unfolding (transforming) bent cells so that their area can be
%calculated with some accuracy

%clear all
%load AD_8c_240to6_60min_19feb_2011_pos_no_2_w_md1

%test cells here are 1 & 4 for tp 140


vol_c =0;%zeros(1,numbM);
area_c=0;%zeros(1,numbM);
min_cell_area=50;

maj_appr_axis=0;%zeros(1,numbM);
min_appr_axis=zeros(1,200);%zeros(200,numbM);

%for t_here=1:numbM

  %   min_appr_axis=0;
%     disp(t_here)
   %no_cells=length(all_obj.max_nucl_int(:,1));
    for j=cell_number %4 %1:no_cells
        
       % tmpI=(all_obj.cells(:,:,t_here)==j);
       tmpI=new_c_Image2; 
        
        if sum(tmpI(:))>min_cell_area % if cell           
            tmpI=bwlabel(tmpI);
         %   index=1;
            %--- remove multiple obj - retain largest
            if max(tmpI(:))>1;
                max_ar=0;index=0;
                for i=1:max(tmpI(:));
                    if sum(sum(tmpI==i))>max_ar;
                        max_ar=sum(sum(tmpI==i)); index=i;
                    end;
                end;
            else
                index=1;
            end;
            tmpI=bwlabel(tmpI==index);
            %--- end removal of multiple objects     
            
            S1=regionprops(bwlabel(tmpI),'boundingbox','area','Orientation','ConvexImage');
            
            xc=round(S1.BoundingBox(1)):-1+round((S1.BoundingBox(1)+S1.BoundingBox(3)));
            yc=round(S1.BoundingBox(2)):-1+round((S1.BoundingBox(2)+S1.BoundingBox(4)));
            tmpI2=tmpI(yc,xc);
            % tmpc1=bwmorph(tmpI2,'remove',inf);
            % tmpc2=bwmorph(tmpI2,'skel',inf);
            % tmpc3=bwmorph(tmpc1,'skel',inf);            
            %magesc(tmpc1+2.*tmpc2);title(j);pause                       
            % imagesc(tmpI2)
            % imagesc(tmpI2+S1.ConvexImage)
            
            B = bwlabel( imrotate(tmpI2,-S1.Orientation));
          %   figure(1);imagesc(B);figure(2);imagesc(tmpI2)           
            % sum(sum(bwmorph(tmpI2,'remove',inf)))
            % sum(sum(bwmorph(B,'remove',inf)))            
           
            %--- remove multiple obj - retain largest         
            if max(B(:))>1;
                max_ar=0;index=0;
                for i=1:max(B(:));
                    if sum(sum(B==i))>max_ar;
                        max_ar=sum(sum(B==i)); index=i;
                    end;
                end;
                  else
                index=1;
            end;
            B=bwlabel(B==index);
            %--- end removal of multiple objects            
            S1=regionprops(bwlabel(B),'boundingbox','area','Orientation','ConvexImage');                    
         %  S1
            xc=round(S1.BoundingBox(1)):-1+round((S1.BoundingBox(1)+S1.BoundingBox(3)));
            yc=round(S1.BoundingBox(2)):-1+round((S1.BoundingBox(2)+S1.BoundingBox(4)));
            mpI2=B(yc,xc);
            [xs,ys]=size(mpI2);
           
            maj_appr_axis=ys;
        %   min_appr_axis=0;
            area_ht=0;
            vol_ht =0;
            for i=1:ys
                r_h=(sum(mpI2(:,i))./2); %factor pi removed by KS
                area_ht=area_ht+2.*pi.*r_h;
                vol_ht=vol_ht+pi.*(r_h).^2;
                
                min_appr_axis(1,i)=sum(mpI2(:,i));
                
            end
         %   vol_c(t_here)=vol_ht;
         %   area_c(t_here)=area_ht;   
            vol_c=vol_ht;
            area_c=area_ht;   
            
        end %if cell for this time
    end %no of cells
end %time

