function shuffled = shuffleSignals(singles)
     [t,n] = size(singles);
     for i=1:n
         shift = round(rand*t);
         temp = singles(:,i);
         singles(:,i) = circshift(temp,shift);
     end
     shuffled = singles;
end