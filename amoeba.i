func amotry(&p,y,&psumm,funcName,ihi,fac)
{
  /* Extrapolates by a factor fac through the face of the simplex, across
     from the high point, tries it and replaces the high point if the new
     point is better.
  */

  fac1 = (1.0 - fac) / numberof(psumm);
  fac2 = fac1  - fac;
  ptry = psumm * fac1 - p(,ihi) * fac2;
  ytry = funcName(ptry); //Eval fcn at trial point
  if (ytry < y(ihi)) {   //If its better than highest, replace highest
    y(ihi) = ytry;
    psumm = psumm + ptry - p(,ihi);
    p(1:,ihi) = ptry;
  }
  return ytry;
}


func amoeba(ftol,funcName,&nCalls,&y,nMax=,p0=,scale=,p=)
{
  if (scale != []) {    //If set, then p0 is initial starting pnt
    ndim = numberof(p0);
    p = p0(,-::ndim);
    for (i=1;i<=ndim;i++) p(i,i+1) = p0(i) + scale(clip(i,,numberof(scale)));
  }

  s = dimsof(p);
  if (s(1) != 2) write, "Either (scale,p0) or p must be initialized";
  ndim = s(2);			//Dimensionality of simplex
  mpts = ndim+1;			//# of points in simplex
  if (nMax == [])  nMax = long(5000);

  val = funcName(p(,1));
  y = array(val, mpts);  //Init Y to proper type
  for (i=2;i<=ndim+1;i++) y(i) = funcName(p(,i));   //Fill in rest of the vals
  nCalls = 0;
  psumm = p(,sum);

  do { //Each iteration
    s = sort(y);
    ilo = s(1);		//Lowest point
    ihi = s(ndim+1);		//Highest point
    inhi = s(ndim);	//Next highest point
    d = abs(y(ihi)) + abs(y(ilo)); //Denominator = interval
    if (d != 0.0) rtol = 2.0 * abs(y(ihi)-y(ilo))/d;
    else rtol = ftol / 2.;         //Terminate if interval is 0

    if (rtol < ftol) {//Done?
      t = y(1);
      y(1) = y(ilo);
      y(ilo) = t;   //Sort so fcn min is 0th elem
      t = p(,ilo);
      p(,ilo) = p(,1);
      p(,1) = t;
      return t;                 //params for fcn min
    }
    
    nCalls = nCalls + 2;
    ytry = amotry(p, y, psumm, funcName, ihi, -1.0);
    if (ytry <= y(ilo)) ytry = amotry(p,y,psumm, funcName,ihi,2.0);
    else if (ytry >= y(inhi)) {
      ysave = y(ihi);
      ytry = amotry(p,y,psumm,funcName, ihi, 0.5);
      if (ytry >= ysave) {
        for (i=1;i<=ndim+1;i++) {
          if (i != ilo) {
            psumm = 0.5 * (p(,i) + p(,ilo));
            p(,i) = psumm;
            y(i) = funcName(psumm);
          }
        }
        nCalls += ndim;
        psumm = p(,sum);
      }		//ytry ge ysave
    } else nCalls -= 1;
  } while (nCalls < nMax);
    
  return -1;		//Here, the function failed to converge.
}

func my_func(p) {
  x=(float(indgen(17))-1.)*5.;
  y=[ 12.0, 24.3, 39.6, 51.0, 66.5, 78.4, 92.7, 107.8, 120.0, 135.5, 147.5, 161.0, 175.4, 187.4, 202.5, 215.4, 229.9];
  return max(abs(y-(p(1)+p(2)*x)));
}

func test_amoeba(void)
// test function similar to the idl one
//IDL prints:
//Intercept, Slope:      11.4100      2.72800 
//Function value:       1.33000 
{
  r=amoeba(1.e-5,my_func,nc,fval,p0=[0.,0.],scale=1.e2);
  r;
  fval;
}
