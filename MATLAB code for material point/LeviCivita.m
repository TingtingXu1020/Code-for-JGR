%Levi - Civita Symbol Implemenetation
%Takes a list of indices [i,j,k,...] and returns the scalar
%Levi-Civita value associated with them, e_i,j,k...
%0 for duplicate entries, otherwise +1/-1 based on ordering or chirality
%
function V = LeviCivita(LArgs)
mp = max(LArgs);
LNot = []; LNC = 0;
%Find Mising Indices between 1 and max index
for I = 1:mp
    if(~isin(I,LArgs))
        LNC = LNC + 1;
    end
end
%Create appended list (if in the other order, then list is skewed)
List = [LNot,LArgs];
%Check for duplicate entries
Dup = 0;
for I = 1:mp
    FoundOnce = 0;
    for J = 1:size(List,2)
        if(and(List(J) == I,FoundOnce == 0))
            FoundOnce = 1;
        elseif(List(J) == I)
            Dup = 1;
        end
        if(Dup == 1)
            break;
        end
    end
    if(Dup == 1)
        break;
    end
end
if(Dup == 1)
    V = 0; %if duplicate entries, LCS = 0
else
   %Count swaps necessary to order the list 
   Count = 0;
   for I = 1:max(List)
       for J = I:size(List,2)
           Breakflag = 0;
           if(List(J) == I)
               Count = Count + J-I; %swaps j into i position
               if((J-I)>1)
                   Count = Count + (J-I)-1; %swaps i into former j position
               end
               List(J) = List(I); %Swap index position
               List(I) = I;
               Breakflag = 1;
           end
           if(Breakflag == 1)
               break;
           end
       end
   end
   %if the swap count is even, LCS = +1, else -1
   if(mod(Count,2)==0)
       V = 1;
   else
       V = -1;
   end
end

return

function A = isin(Val,List)
A = 0;
for I = 1:size(List,2)
    if(List(I) == Val)
        A = 1;
    end
end
return;
