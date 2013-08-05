void shell(unsigned long n, double a[])
     /* sorts an array a[1..n] into ascending numerical order by Shell's
	method */
{
  unsigned long i,k,inc;
  double v;

  inc=1;
  do{
    inc*=3;
    inc++;
  }while(inc<=n);
  do{      /* loop over partial sorts */
    inc/=3;
    for(i=inc+1;i<=n;i++){
      v=a[i];
      k=i;
      while(a[k-inc]>v){  /* inner loop */
	a[k]=a[k-inc];
	k-=inc;
	if(k<=inc) break;
      }
      a[k]=v;
    }
  }while(inc>1);
}
