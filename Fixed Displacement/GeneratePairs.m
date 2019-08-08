items=[0.56,0.23,0.23, 0.23, 0.42, 0.62, 0.62,1.11, 0.70,0.33,0.52,0.52,1.52,...
    0.52,0.23,0.23,0.46,0.96,0.23,1.10, 1.8, 1.71, 0.87, 1.34];

dblUniqItems=[unique(items),unique(items)];
items= sort(items);

fits=[3.73,3.76];
% 
% v=findCombination(1.9,items)
% sum(v)

% all on one side
% 3.74 fits
% 3.75 barely fits
% 3.76 barely fits
% 3.77 doesnt really fit
% 3.73 fits
% 3.72, cant be created
% 3.71, cant be done
% 3.70 too loose.
% 3.8 doesnt fit

% items= sort(dblUniqItems);


totalDisplacement=@(x) 3.76+x*0.0737;
TargetDisplacements=0:0.1:1.9;
p1Save=zeros(20,5);
p2Save=zeros(20,5);
for n =1:size(TargetDisplacements,2)
%    find first combination for one side
    p1=findCombination(TargetDisplacements(n),items,true);
%     create a list that excludes the used values
   tempItems=items;
   for m=1:length(p1)
        tempItems(find(tempItems==p1(m),1))=[];
   end
%    find other side
   p2=findCombination(totalDisplacement(TargetDisplacements(n))-sum(p1),tempItems,false);
   
   p1Error(n)=TargetDisplacements(n)-sum(p1);
   TotalError(n)=totalDisplacement(TargetDisplacements(n))-(sum(p1)+sum(p2));
%    if abs(p1Error(n))>0.03 | TotalError(n)>0.03
%        p1
%        p2
%    end
    for m=1:length(p1)
        p1Save(n,m)=p1(m);
    end
    for m=1:length(p2)
        p2Save(n,m)=p2(m);
    end
   p1;
   p2;
end

for n=size(TargetDisplacements,2):-1:1
    p1Save(n,:)
    p2Save(n,:)
end