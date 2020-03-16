clear; close all; clc; 


promp = "Enter x quantity:  ";
xQuantity = input(promp)
c = zeros(2*xQuantity);
c=c(1,:);
A = zeros(xQuantity, xQuantity);
b = zeros(xQuantity);
b=b(:,1);
I = eye(xQuantity);

for i=1:xQuantity 
    promp = "Enter coefficiant for x" + i + ": ";
    c(i) = input(promp); 
end

promp = "Press: 1 if MAX ----- 2 if MIN  ";
choice = input(promp);
sign = zeros(xQuantity);

for i=1:xQuantity
   promp = i +". border condition:"
   for j = 1:xQuantity
        promp2 = "__ x" + j + " --> ";
        A(i,j) = input(promp2);
   end
            promp3 = "Press: 1 if <= ----- 2 if >=  ";
            sign(i) = input(promp3); 
            disp("---------------")
            promp2 = "Enter second side of equation ";
            b(i) = input(promp2);
end


for i=1:xQuantity
   if (choice == 1)
       if (sign(i)==2)
           b(i) = -b(i);
           A(i,:) = - A(:,i);
       end
   else
       if (sign(i)== 1)
           b(i) = -b(i);
           A(:,i) = - A(:,i);
       end 
   end
end

nextA = A; 
nextb = b;
cb = zeros(xQuantity);
cb=cb(1,:);
xb = zeros(xQuantity);
xb=xb(1,:);

if(choice == 1)
    for i=1:xQuantity
        xb(i)=xQuantity+i; 
    end
else
    for i=1:xQuantity
        xb(i)=2*xQuantity+i; 
    end
end 

B = I;

disp("****************** CANONICAL FORM: *******************")
if (choice == 1)
    z = zeros(2*xQuantity);
    z=z(1,:);
    simplexCriterion = zeros(2*xQuantity);
    simplexCriterion = simplexCriterion(1,:);
     disp("Target Function:");
        for i=1:2 * xQuantity
            if (i <= xQuantity) 
                fprintf(c(i) + " x" + i + " + ");
            elseif (i == 2 * xQuantity) 
                    fprintf("0 x" + i + " --> MAX ");
            else
                 fprintf("0 x" + i + " + ");
            end
        end
           
       
        disp("");
        disp("Border Conditions: ");
        for i=1:xQuantity
            disp("  ");
            disp(i + ". BC :  ");
            for j=1:xQuantity
               fprintf( A(i,j) + " x" + j + " + ");
                if (j == xQuantity) 
                    fprintf(1.0 + " x" + (i + xQuantity) + " = " + b(i));
                end
            end
            disp("");
        end
        disp("  ");
        
        stop = false;
        iteration = 1;
        
    while (~stop)
        disp(" ************ Simplex Matrix No. " + iteration + " ************") 
       invB = inv(B); 
       joined = [nextA invB];
       z = cb*joined;
       
        for i=1:2*xQuantity
            simplexCriterion(i) = c(i) - z(i);
        end
        
                
    
    [maxValue, inside] = max(simplexCriterion);
    divided =b./joined(:,inside);
    [minValue, outside] = min(divided);      
    outside;
    
    disp("cb")
    cb'
    disp("xb --> x in base")
    xb'
    joined
    nextb
    divided
    %X=[cb' xb' joined nextb divided]
    z
    simplexCriterion
    
        more = 0;
        for i=1:2*xQuantity 
                if (simplexCriterion(i) > 0) 
                    more = more+1;  
                end
        end
         
        if (more <= 0 ) 
            stop = true;
        end
        
        
            disp("_____________________________");
            disp("X" + xb(outside) + " goes out");
            disp("X" + (inside) + " goes in");

            cb(outside) = c(inside);
            xb(outside) = inside;

            for i=1:xQuantity
                B(i,outside) = A(i,inside);
            end
        
        nextA = inv(B)*A;
        nextb = inv(B)*b;
        iteration=iteration+1;       
                
    end
    
    
   solution = sum(cb * nextb) 
        
        
else
    z = zeros(3*xQuantity);
    z=z(1,:);
    simplexCriterion = zeros(3*xQuantity);
    simplexCriterion = simplexCriterion(1,:);
    
    for i=1:xQuantity 
        c(i+2*xQuantity) = 1000; 
        % can't choose Inf, computer can't calculate for example (0*Inf) 
        cb(i) = 1000; 
    end
    
    
    disp("Target Function:");
        for i=1:3 * xQuantity
            if (i <= xQuantity) 
                fprintf(c(i) + " x" + i + " + ");
            elseif (i <= 2*xQuantity)
                 fprintf("0 x" + i + " + ");
            elseif (i <3*xQuantity)
                fprintf("m x" + i + " + ");
            else
                fprintf("m x" + i + " --> MIN ");
            
            end
        end
        disp("")
        
        disp("Border Conditions: ");
        for i=1:xQuantity
            disp("  ");
            disp(i + ". BC :  ");
            for j=1:xQuantity
               fprintf( A(i,j) + " x" + j + " + ");
                if (j == xQuantity) 
                    fprintf("(" + -1.0 + ")" + " x" + (i + xQuantity) + " + " );
                    fprintf(1.0 + " x" + (i + 2*xQuantity)  + " = " + b(i));
                end
            end
            disp("");
        end
        disp("  ");
        stop = false;
        iteration = 1;
        helpB = -eye(xQuantity)
        tmp = 0; 
        inside = -10;
        
    while (~stop && tmp ~= 2 )
        
        disp(" ************ Simplex Matrix No. " + iteration + " ************") 
       invB = inv(B); 
       joined = [nextA -invB invB];
       z = cb*joined;
       
       
        for i=1:3*xQuantity
            simplexCriterion(i) = c(i) - z(i);
        end
        
                
    prev_inside = inside; 
    [minValue2, inside] = min(simplexCriterion);
    
    
    divided = nextb./joined(:,inside);
    tmp2 = divided; 
    tmp2(divided<=0) = NaN  %Set anything <= 0 to NaN
    [minValue, outside] = min(tmp2);      
    
    if (inside == xb(outside) | prev_inside == xb(outside)) 
          tmp = tmp +1;
    end
    
    disp("cb")
    cb'
    disp("xb --> x in base")
    xb'
    joined
    nextb
    divided
    z
    simplexCriterion
    
        more = 0;
        for i=1:3*xQuantity 
                if (simplexCriterion(i) < 0) 
                    more = more+1;  
                end
        end
         
        if (more <= 0 ) 
            stop = true;
        end
        
        
            
            disp("_____________________________");
            disp("X" + xb(outside) + " goes out");
            disp("X" + (inside) + " goes in");
            
            
            cb(outside) = c(inside);
            xb(outside) = inside;
            
            
            
           

            for i=1:xQuantity
                if(inside<= xQuantity)
                    B(i,outside) = A(i,inside);
                else
                    B(i,outside) = helpB(i,inside-xQuantity);
                end
            end
            
        
        nextA = inv(B)*A;
        nextb = inv(B)*b;
        iteration=iteration+1;  
        
                
    end
    
    
   solution = sum(cb * nextb) 
   
   
        
end





