function results=run_all(D)
if nargin<1,D=24;end
results=cell(3,2,4);
lengths=[128];kinds={'polar','pac'};modes={'none'};
for i=1:1
    for j=1:2
        for k=1:1,results{i,j,k}=run_spectrum(lengths(i),kinds{j},modes{k},D);end
    end
end
end
