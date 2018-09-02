function cnt = stanford_graphs( n,graph,pos_out )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
[G,delout] = importdata(graph);

for j = 1:length(G)
   G(j,1) = G(j,1)+1;
   G(j,2) = G(j,2)+1;
end

for j = 1:n
   L{j} = [];
end

for j = 1:length(G)
   L{G(j,1)} = [L{G(j,1)} G(j,2)];

end

for j = 1:n
   c(j) = length(L{j});
end

% Power method

p = .85;
delta = (1-p)/n;
x = ones(n,1)/n;
z = zeros(n,1);
cnt = 0;
tic
while max(abs(x-z)) > .0000001
   z = x;
   x = zeros(n,1);
   for j = 1:n
      if c(j) == 0
         x = x + z(j)/n;
      else
         x(L{j}) = x(L{j}) + z(j)/c(j);
      end
   end
   x = p*x + delta;
   cnt = cnt+1;
end
toc


dlmwrite(pos_out,x,'delimiter',' ');

end

