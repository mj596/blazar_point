int 
locate(double xx[], unsigned long n, double x)
     /* Given an array xx[1..n], and given a value x, returns a value i
	such that x is between xx[i] and xx[i+1]. xx must be increasing.
	i=0 or i=n is returned to indicate that x is out of range. */
{
  unsigned long iu,im,il;
  
  il=0;
  iu=n+1;
  while(iu-il > 1){
    im=(iu+il) >> 1;
    if(x>=xx[im])
      il=im;
    else
      iu=im;
  }
  return il;
}
