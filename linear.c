double linear(x,x1,y1,x2,y2)
double x,x1,y1,x2,y2;
{
    double a,b;

    a=(y1-y2)/(x1-x2);
    b=0.5*(y1+y2-(y1-y2)*(x1+x2)/(x1-x2));
    return(a*x+b);
}
