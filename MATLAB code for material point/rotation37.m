function[inputrotation]=rotation37(numc)
% catetian coordinates
%Bazant1986
%x1=cos(phi)sin(theta); x2=sin(phi)sin(theta); x3=cos(theta)
inputorien=zeros(3, numc);

% spherical coordinates
inputspherical=ones(3, numc);

% rotation of the orientation vector
inputrotation=zeros(3, numc);

inputorien(:,1)=[1 0 0]';
inputorien(:,2)=[0 1 0]';
inputorien(:,3)=[0 0 1]';
inputorien(:,4)=[0.707106781 0.707106781 0]';
inputorien(:,5)=[0.707106781 -0.707106781 0]';
inputorien(:,6)=[0.707106781 0 0.707106781]';
inputorien(:,7)=[0.707106781 0 -0.707106781]';
inputorien(:,8)=[0 0.707106781 0.707106781]';
inputorien(:,9)=[0 0.707106781 -0.707106781]';
inputorien(:,10)=[0.951077869651 0.308951267775 0]';
inputorien(:,11)=[0.951077869651 -0.308951267775 0]';
inputorien(:,12)=[0.308951267775 0.951077869651 0]';
inputorien(:,13)=[0.308951267775 -0.951077869651 0]';
inputorien(:,14)=[0.951077869651 0 0.308951267775]';
inputorien(:,15)=[0.951077869651 0 -0.308951267775]';
inputorien(:,16)=[0.308951267775 0 0.951077869651]';
inputorien(:,17)=[0.308951267775 0 -0.951077869651]';
inputorien(:,18)=[0 0.951077869651 0.308951267775]';
inputorien(:,19)=[0 0.951077869651 -0.308951267775]';
inputorien(:,20)=[0 0.308951267775 0.951077869651]';
inputorien(:,21)=[0 0.308951267775 -0.951077869651]';
inputorien(:,22)=[0.335154591939 0.335154591939 0.880535518310]';
inputorien(:,23)=[0.335154591939 0.335154591939 -0.880535518310]';
inputorien(:,24)=[0.335154591939 -0.335154591939 0.880535518310]';
inputorien(:,25)=[0.335154591939 -0.335154591939 -0.880535518310]';
inputorien(:,26)=[0.335154591939 0.880535518310 0.335154591939]';
inputorien(:,27)=[0.335154591939 0.880535518310 -0.335154591939]';
inputorien(:,28)=[0.335154591939 -0.880535518310 0.335154591939]';
inputorien(:,29)=[0.335154591939 -0.880535518310 -0.335154591939]';
inputorien(:,30)=[0.880535518310 0.335154591939 0.335154591939]';
inputorien(:,31)=[0.880535518310 0.335154591939 -0.335154591939]';
inputorien(:,32)=[0.880535518310 -0.335154591939 0.335154591939]';
inputorien(:,33)=[0.880535518310 -0.335154591939 -0.335154591939]';
inputorien(:,34)=[0.577350269190 0.577350269190  0.577350269190]';
inputorien(:,35)=[0.577350269190 0.577350269190  -0.577350269190]';
inputorien(:,36)=[0.577350269190 -0.577350269190  0.577350269190]';
inputorien(:,37)=[0.577350269190 -0.577350269190  -0.577350269190]';

for i=1:numc
    inputspherical(2, i) =acos(inputorien(3,i)); %theta
    if inputspherical(2, i) == 0 ||  inputorien(1, i)==0 % theta = 0
        inputspherical(3, i) = 0;
    else
       inputspherical(3, i) = atan(inputorien(2,i)/inputorien(1,i));
    end
%     inputspherical(3, i)=atan(inputorien(2,i)/inputorien(1,i)); %phi
%      if inputorien(1, i)==0 && inputorien(2,i)==0
%         inputspherical(3, i)=0;%phi
%      end
    inputrotation(1,i)=inputspherical(3,i); %phi
    inputrotation(2,i)=inputspherical(2,i);  %theta
end


end