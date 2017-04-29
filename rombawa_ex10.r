Mullers = function(f,x0,x1,x2,macheps,max,verbose) #function, first three estimates, machine epsilon, max iterations, verbose
{
  #0. set up default values
  if(missing(max))
  {
    max = 100000;
  }
  
  if(missing(verbose))
  {
    verbose = TRUE;
  }
  
  if(missing(macheps))
  {
    macheps = 0.000000001; #1*10^-9 unless naduling ako hahaha
  }
  count = 1;
  x3 = 0;
  
  h0 = 0;
  h1 = 0;
  d0 = 0;
  d1 = 0;
  
  a = 0;
  b = 0;
  c = 0;
  
  checkpos = 0; #holder for |b+sqrt(b^2-4ac)|
  checkneg = 0; #holder for |b-sqrt(b^2-4ac)|
  sign = FALSE; #true if positive, false if negative
  
  #iterative computation
  
  while(count <= max)
  {
    #1. compute for h0,h1,d0,d1
    
    fx0 = f(x0);
    fx1 = f(x1);
    fx2 = f(x2);
    
    h0 = x1-x0;
    h1 = x2-x1;
    d0 = ( f(x1)-f(x0) ) / (x1-x0);
    d1 = ( f(x2)-f(x1) ) / (x2-x1);
    
    #2. compute for a,b,c
    a = (d1-d0) / (h1-h0);
    b = a*h1 + d1;
    c = fx2;
    
    #3. check for sign
    checkpos = b + sqrt( b*b - 4*a*c );
    checkneg = b - sqrt( b*b - 4*a*c );
    
    if(abs(checkneg) < abs(checkpos))
    {
      sign = TRUE; #positive
    }
    else
    {
      sign = FALSE;
    }
    
    #4. compute for x3, using form according to sign
    if(sign == TRUE)
    {
      x3 = x2 - (2*c) / checkpos;
    }
    else
    {
      x3 = x2 - (2*c) / checkneg;
    }
    
    #5. compute approximate error
    ea = abs((x3-x2)/x3);
    
    #6. if verbose, print values
    if(verbose == TRUE)
    {
      # print: count, x0, x1, x2, x3, ea
      print("Iteration:");
      print(count);
      print("x0");
      print(x0);
      print("x1");
      print(x1);
      print("x2");
      print(x2);
      print("Root Estimate: x3");
      print(x3);
      print("Approximate Relative Error:");
      print(ea);
    }
    
    #7. terminating condition: error less than machine epsilon
    if(ea <= macheps)
    {
      return(x3);
    }
    
    #8. adjust values for next iteration
    x0 = x1;
    x1 = x2;
    x2 = x3;
    
    count = count+1;
  }
  #worst case scenario: count maxed out
  print("Maximum iterations exhausted");
  return(x3);
}

# test scenario

f = function(x) return(x*x - 4);
x0 = 0;
x1 = 1;
x2 = 4;
max = 100000;
verbose = FALSE;
macheps = 0.00001;
solution = Mullers(f,x0,x1,x2);
print("SOLUTION")
print(solution);