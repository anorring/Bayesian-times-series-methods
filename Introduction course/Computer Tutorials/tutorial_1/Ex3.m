    %This is a program which illustrates do loops and if statements 
    %first create a column vector to work with
    x=[1;2;7;5;9;3;6;9;1;11;1];
    %the following command sums up the elements of a column vector
    xsum=sum(x);
    xsum1=0;
    for i=1:11
        xsum1=xsum1 + x(i,1);
    end
    disp xsum;
    disp(xsum);
    disp xsum1;
    disp(xsum1);
    %now illustrate the if command
        xsum2=0;
    for i=1:11
    if x(i,1)>4
    xsum2=xsum2 + x(i,1);
    end
    end
    disp xsum2
    disp(xsum2);
    
