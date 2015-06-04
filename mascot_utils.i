func does_it_fit(slist)
{
  nstars = dimsof(slist)(0);
  if (nstars==1) return 1;

  if (nstars==2) {
    d = slist(1:2,1)-slist(1:2,2);
    d = sqrt(sum(d^2.));
    if (d<=(120-2*edge_margin)) return 1;
    else return 0;
  }

  if (nstars==3) {
    d = array(0.,[3,300,300,nstars]);
    for (ns=1;ns<=3;ns++) {
      d(,,ns) = dist(300,xc=150+slist(1,ns),yc=150+slist(2,ns));
    }
    dmin = min(d(,,max));
    if (debug>2) write,format="dmin=%f\n",dmin;
    if (dmin<=(60-edge_margin)) return 1;
    else return 0;
  }
  return 0; // shouldn't be here anyway.
}

func escapechar(s,eseq)
{
  if (eseq==[]) eseq="!"; // primarily for yorick graphics
  s=streplace(s,strfind("_",s,n=20),eseq+"_");
  s=streplace(s,strfind("^",s,n=20),eseq+"^");
  return s;
}


func wfs_noise(mag,verbose=)
/* DOCUMENT wfs_noise(mag,verbose=)
   Return the noise on the TT sensors in arcsec rms.
   SEE ALSO:
 */
{
  nstars = numberof(mag);
  ttnpde  = array(0.,nstars);
  for (i=1;i<=numberof(mag);i++) ttnpde(i) = magstar(mag(i));

  ttnpde *= detqe*thrup/sampfreq;
  ttnsky  = (magsky())(case_sky)*detqe*thrup/sampfreq/4.; // 4-> per channel
  ttdc    = dark_current/sampfreq;
  ttsnr   = ttnpde/sqrt(ttnpde+4*ttron^2.+4*ttnsky+4.*ttdc);
  if (correct_for_poisson) {
    ttsnr = ttsnr/(1.+(0.7+0.4*sqrt(ttnsky+ttdc))/ttnpde)
  }
  // FIX ME ! don't know what gain FwhmCloseLoop is ...
  // gainfwhmcloseloop = 1.8; // how much the FWHM is improved at lambdawfs in CL
  gainfwhmcloseloop = 1.; // let's be conservative and assume no gain !
  r0wfs = r0vis * (lambdawfs/0.5)^1.2;
  
  // rms TT noise in as on WFSs:
  // the 0.533 is for the CG (CG=FWHM*0.533)
  ttnoise= 0.533*(lambdawfs*1.0e-6/r0wfs/gainfwhmcloseloop/4.848e-6)/ttsnr;

  if (verbose!=0) {
    write,format="%s","TT npde[ph/frame] = "; write,ttnpde;
    write,format="TT Nsky[ph/channel/frame] = %.2f  Dark/channel/frame = %.2f\n",ttnsky,ttdc;
    write,format="%s","TT noise[as] = "; write,ttnoise;
  }

  ttnv = [];
  for (i=1;i<=nstars;i++) grow,ttnv,_(ttnoise(i),ttnoise(i));
  return ttnv; // 2 * nstars vector, sqrt(diagonal of noise covar mat).
}


func ftcb(te,tcal,tmir,gain,dim,x=)
/* DOCUMENT ftcb(te,tcal,tmir,gain,dim,x=)
   returns [f,hbo,hcor,hbf]
   AUTHOR: F.Rigaut, way back in 1996?
   SEE ALSO: 
 */
{
  f = indgen(dim)/te/2./dim;
  if (!is_void(x)) { f = x;}
  p = 2i*pi*f;

  hzoh  = (1.-exp(-te*p))/(te*p);
  hmir  = 1./(1.+tmir*p);
  hwfs  = (1.-exp(-te*p))/(te*p);
  hcal  = gain*exp(-tcal*p);

  hbo  = hzoh*hmir*hwfs*hcal/(1-exp(-p*te));

  hcor = double((1./(1.+hbo))*conj(1./(1.+hbo)));
  hbf  = double((hbo/(1.+hbo))*conj(hbo/(1.+hbo)));
  hbo  = double(hbo*conj(hbo));

  return ([f,hbo,hcor,hbf]);
}



func plvf(vy,vx,y,x,autoscale=,scale=,width=,hsize=,hang=,color=,type=,prop=)
/* DOCUMENT plvf,vy,vx,y,x,scale=,width=,hsize=,hang=,color=,type=,prop=
   Plots the vector field defined by (vx,vy) at positions (x,y)
   vx,vy,x,y must have the same size, but can be of arbitrary dimension.
   KEYWORDS:
   autoscale: set to 1 for the vector length to be autoscaled
   scale:     multiplicative factor applied to the autoscale results
              (for fine tweaking)
   width, color, type: same as in plg.
   hsize, hang: size and opening angle of the arrow head (default
       hsize=0.4, hang=20 degrees)
   prop:      set to zero if you want the same head size for *all* vector.
              Otherwise, the head size is proportionnal to the size of
              the vector (which results in something nicer to the eye).
   SEE ALSO: pldj
 */
{
  vx = vx(*); vy = vy(*); x = x(*); y = y(*);
  
  if (!scale) scale=1.;
  if (!width) width=2;
  if (hsize==[]) hsize=0.4;
  if (hang==[]) hang = 20;
  if (prop==[]) prop = 1;

  if (autoscale) {  
    if (prop) {
      sc=abs(vx,vy);
      if (max(sc)==0) sc=1.;
      //      else sc=sc/max(sc);
    } else {sc=1.;}

    // vector body autoscaling:
    xdif = abs(x(dif));
    w = where(xdif != 0);
    if (numberof(w)!=0) {
      minspace = min(xdif(w));
    }
    ydif = abs(y(dif));
    w = where(ydif != 0);
    if (numberof(w)!=0) {
      minspace = (minspace==[]? min(ydif(w)) : min([minspace,min(ydif(w))]) );
    }
    if (minspace==[]) minspace=1.;
    // autoscale normalization factor: max vector length / min space between location
    norm = max(abs([vy,vx]))/minspace*1.2;
    if (norm==0) norm=1.;
    vx = vx/norm*scale;
    vy = vy/norm*scale;
    //    hsize = hsize/norm*scale;
  } else {
  }
  sc = abs(vx,vy);

  pldj,(x+vx)(*),(y+vy)(*),x(*),y(*),width=width,color=color,type=type;
  x1=(x+vx)(*);  y1=(y+vy)(*);
  ang=atan(vy(*),vx(*))-(180-hang)*pi/180.;
  x2=x1+sc*hsize*cos(ang);
  y2=y1+sc*hsize*sin(ang);
  pldj,x2,y2,x1,y1,width=width,color=color,type=type;

  ang=atan(vy,vx)-(180+hang)*pi/180.;
  x2=x1+sc*hsize*cos(ang);
  y2=y1+sc*hsize*sin(ang);
  pldj,x2,y2,x1,y1,width=width,color=color,type=type;
}




func zernike_spectra(nzernike,&x,&y)
{
  extern nu,v,nzer,deg,rad,tcos;

  require,"qromo.i";

  nzer = nzernike;
  nm   = zernumero(nzer);
  deg  = nm(1);
  rad  = nm(2);
  tcos = -1;

  if (rad != 0) {
    if ((nzer % 2) == 0) tcos=1;
    else tcos=0;
  }

  write,"Zernike number : ",nzer;
  write,"Radial order   : ",rad;
  if (tcos == 1) write,"Cosinus term";
  if (tcos == 0) write,"Sinus term";

  nu    = 1.;
  v     = 1.;

  x     = 0.;
  y     = 0.;

  for (nu=0.001;nu<=20;nu+=0.01) {
    res = qromo(myfunc,0.,20.)*2.;
    //write,res;
    grow,x,nu;
    grow,y,res;
  }

  x     = x(2:);
  y     = y(2:);
}

func myfunc(fy)
{
  extern nu,v,nzer,deg,rad,tcos;

  fx    = nu/v;
  theta = atan(fy,fx);
  sincos = tcos*cos(rad*theta)+(1.-tcos)*sin(rad*theta);
  if (rad == 0) sincos = 1;

  val   = (1./v)*(fx^2.+fy^2.)^(-17./6.)*\
    (abs(bessj(deg+1,2*pi*sqrt(fx^2.+fy^2.))))^2.*sincos^2.;

  return val;
}

func zernumero(zn)
{
  j=0;
  for (n=0;n<=100;n++) {
    for (m=0;m<=n;m++) {
      if ((n-m) % 2 == 0) {
        j=j+1;
        if (j == zn) return [n,m];
        if (m != 0) {
          j=j+1;
          if (j == zn) return [n,m];
        }
      }
    }
  }
}



func vib_spectra(void)
{
  spv = array(0.,[2,4000,6]);
  spv(,1) = freq = span(0.,1000.,4000);

  for (i=1;i<=numberof(*tipvibfreq);i++) {
    spv(,2) += ((*tipvibrms)(i))*                                       \
      exp(-((freq-(*tipvibfreq)(i))/((*tipvibwidth)(i)/1.66))^2.)^2.;
  }
  spv(,2) = spv(,2)/sum(spv(,2));

  for (i=1;i<=numberof(*tiltvibfreq);i++) {
    spv(,3) += ((*tiltvibrms)(i))*                                       \
      exp(-((freq-(*tiltvibfreq)(i))/((*tiltvibwidth)(i)/1.66))^2.)^2.;
  }
  spv(,3) = spv(,3)/sum(spv(,3));

  return spv;
}
  
func null_modes_spectra(void)
{
  sp = array(0.,[2,4000,6]);

  if (fileExist("zernike_spectra.fits")) {
      write,"Reading zernike_spectra.fits";
      a = fits_read("zernike_spectra.fits");
      sp(1:numberof(a(,1)),) = a;
  } else {
    zernike_spectra,2,x2,y2;
    a=array(0.,[2,numberof(x2),6]);
    a(,1)=x2; a(,2)=y2;
    zernike_spectra,3,x3,y3;
    a(,3)=y3;
    zernike_spectra,4,x4,y4;
    a(,4)=y4;
    zernike_spectra,5,x5,y5;
    a(,5)=y5;
    zernike_spectra,6,x6,y6;
    a(,6)=y6;
    fitsWrite,"zernike_spectra.fits",a;
    sp(1:numberof(x2),) = fits_read("zernike_spectra.fits");
  }

  //read out zernike spectra
  sp(1:,1)      = double(indgen(4000)-1.)*(sp(3,1)-sp(2,1))+sp(1,1);
  //fill in frequency vector
  freq          = sp(,1);
  dfreq         = (sp(3,1)-sp(2,1));

  rtel          = tel_diam/2.;
  dr0           = tel_diam/(r0vis*(lambdawfs/0.5)^1.2);
  cn2           = cn2/sum(cn2);
  dr0i          = (cn2*dr0^(5./3.))^(3./5.);
  nlayers       = numberof(cn2);

  b = array(0.,[2,4000,6]);
  b(,1) = freq;

  //Tip:
  sp2           = sp(,1)*0.;
  x2 = a(,1);
  dfreqinit = x2(3)-x2(2);
  y2 = a(,2)/(sum(a(,2))*dfreqinit);
  for (i=1;i<=nlayers;i++) {
    tmp         = spline(x2*y2,x2,freq/(wind(i)/rtel))/freq;
    tmp         = tmp/(sum(tmp)*dfreq);
    sp2         = sp2+tmp*0.45*dr0i(i)^(5./3.);
  }
  b(,2) = sp2;

  //Tilt:
  sp3           = sp(,1)*0.;
  y3 = a(,3)/(sum(a(,3))*dfreqinit);
  for (i=1;i<=nlayers;i++) {
    tmp         = spline(x2*y3,x2,freq/(wind(i)/rtel))/freq;
    tmp         = tmp/(sum(tmp)*dfreq);
    sp3         = sp3+tmp*0.45*dr0i(i)^(5./3.);
  }
  b(,3) = sp3;

  //focus:
  sp4           = sp(,1)*0.;
  y4 = a(,4)/(sum(a(,4))*dfreqinit);
  for (i=2;i<=nlayers;i++) {
    tmp         = spline(x2*y4,x2,freq/(wind(i)/rtel))/freq;
    tmp         = tmp/(sum(tmp)*dfreq);
    sp4         = sp4+tmp*alt(i)^2.*0.02332*dr0i(i)^(5./3.);
  }
  b(,4) = sp4;

  //astig:
  sp5           = sp(,1)*0.;
  y5 = a(,5)/(sum(a(,5))*dfreqinit);
  for (i=2;i<=nlayers;i++) {
    tmp         = spline(x2*y5,x2,freq/(wind(i)/rtel))/freq;
    tmp         = tmp/(sum(tmp)*dfreq);
    sp5         = sp5+tmp*alt(i)^2.*0.02332*dr0i(i)^(5./3.);
  }
  b(,5) = sp5;

  //astig:
  sp6           = sp(,1)*0.;
  y6 = a(,6)/(sum(a(,6))*dfreqinit);
  for (i=2;i<=nlayers;i++) {
    tmp         = spline(x2*y6,x2,freq/(wind(i)/rtel))/freq;
    tmp         = tmp/(sum(tmp)*dfreq);
    sp6         = sp6+tmp*alt(i)^2.*0.02332*dr0i(i)^(5./3.);
  }
  b(,6) = sp6;

  sp = b;
  b  = [];
  
  ttsp = sp(,2:3)(,avg); ttsp /= sum(ttsp);
  tasp = sp(,4:6)(,avg); tasp /= sum(tasp);
  sp(,2:3) = ttsp(,-);
  sp(,4:6) = tasp(,-);

  return sp;
}


func magstar(vstar)
{
  lambda = [3.5,4, 5, 6, 7, 8, 9,10,11]*1e-7;      // meters
  qe     = [0., 5,30,37,37,32,23,10,0.]*0.01*1.6;  // percent
  h      = 6.62e-34;
  c      = 3e8;

  // for band u,b,v,r,i
  cw     = [365,440,550,700,900]*1e-9;             // central wavelength
  dw     = [68,98,89,22,24]*1e-3;                  // delta_lambda in microns
  zp     = [-11.37,-11.18,-11.42,-11.76,-12.08];
  // photometric zeropoints in W/cm2/mic
  zpv    = vstar;                                  // Vmag of star
  msky   = zpv+[-0.63,-0.58,0.,-0.52,-0.93];       // G0 star magnitudes

  // m = -2.5 (log10(f) +zp)

  lvec   = double(indgen(10)-1)/9.*750.e-9+350.e-9; // lambda vector
  qeapd  = spline(qe,lambda,lvec);

  f      = 10.^(-0.4*msky+zp);                     // f en W/cm2/mic
  f      = f/(h*c/cw)*pi*(400.^2-50.^2);
  // f en N_photon/pupille_gemini/s/mic
  f      = f*0.1;
  // f en N_photon/pupille_gemini/s/100nm
  tab1   = f; grow,tab1,f(5);
  tab2   = cw; grow,tab2,1100e-9;
  sp     = tspline(10,tab1,tab2,lvec);

  return zero_point_fudge*sum(sp*qeapd);
}


func magsky(void)
{
  lambda   = [3.5,4, 5, 6, 7, 8, 9,10,11]*1e-7;     // meters
  qe       = [0., 5,30,37,37,32,23,10,0.]*0.01*1.6; // percent
  h        = 6.62e-34;
  c        = 3e8;

  // for band u,b,v,r,i
  cw       = [365,440,550,700,900]*1e-9;            // central wavelength
  dw       = [68,98,89,22,24]*1e-3;                 // delta_lambda in microns
  zp       = [-11.37,-11.18,-11.42,-11.76,-12.08];
  // photometric zeropoints in W/cm2/mic
  zpv      = 21.8;                                  // mag V per square arcsec at darkest
  msky     = array(double,[2,5,4]);
  msky(,1) = zpv+[0.  ,0.8 ,0.,-0.9,-1.9];          // darkest
  msky(,2) = zpv+[-1.5,0.2 ,0.,-0.8,-1.6]-0.6;      // 50%
  msky(,3) = zpv+[-2.2,0.  ,0.,-0.4,-0.8]-1.8;      // 80%
  msky(,4) = zpv+[-3. ,-0.5,0.,-0.1,-0.2]-3.3;      // bright

  // m = -2.5 (log10(f) +zp)

  lvec  = double(indgen(10)-1)/9.*750e-9+350e-9;     // lambda vector
  qeapd = spline(qe,lambda,lvec);

  res   = 0.;

  for (i=1;i<=4;i++) {
    f   = 10.^(-0.4*msky(,i)+zp);                   // f in W/cm2/mic
    f   = f/(h*c/cw)*pi*(400.^2-50.^2);             // f in N_photon/pup_gemini/s/mic
    f   = f*0.1;                                    //f in N_photon/pup_gemini/s/100nm
    tab1=f; tab2=cw;
    grow,tab1,f(5); grow,tab2,1100e-9;
    sp  = tspline(8,tab1,tab2,lvec);
    grow,res,sum(sp*qeapd);
  }

  res *= (pi*ttwfs_aper_radius^2.); // for the TT wfs field stop
  
  return zero_point_fudge*res(2:);
}


func ra_hms2decimal(hh,mm,ss)
{
  if (mm==[]) {
    mm = hh(..,2);
    ss = hh(..,3);
    hh = hh(..,1);
  }
  hh = hh+mm/60.+ss/3600.;
  return hh*15.;
}

func ra_decimal2hms(ra)
{
  local hh,mm,ss;
  // while (ra<0) ra += 360.
  hh = ra/15.;
  mm = (hh-long(hh))*60.;
  ss = (mm-long(mm))*60.;
  hh = long(hh); mm = long(mm);

  return _(hh*1.,mm*1.,ss);
}

func dec_dms2decimal(dd,mm,ss)
{
  if (mm==[]) {
    mm = dd(..,2);
    ss = dd(..,3);
    dd = dd(..,1);    
  }
  dds = sign(dd);
  dd = abs(dd)+mm/60.+ss/3600.;
  return dds*dd;
}

func dec_decimal2dms(dec)
{
  local dd,mm,ss;
  dds = sign(dec);
  dd = abs(dec);
  mm = (dd-long(dd))*60.;
  ss = (mm-long(mm))*60.;
  dd = long(dd); mm = long(mm);

  return _(dds*dd*1.,mm*1.,ss);
}




func print_asterism(n,file)
/* DOCUMENT print_asterism(n,file)
   Print out a group of guide stars, with coordinates, magnitude and
   distance to field center.
   Print out either in terminal (file == []) or in a file (file=name
   or =1, in which case a file name in the form of the session
   timestamp will be used.
   SEE ALSO:
 */
{
  local f;
  extern current_objname;
  extern session_filename;
  
  nast = numberof(sall);

  if (n==[]) n = current_ast;
  
  if (n>nast) exit,"Asterism# > available number";
  if (n<0)    exit,"Invalid asterism";

  if (file!=[]) {
    if (structof(file)!=string) {
      if (session_filename==[]) {
        tmp = timestamp();
        tmp = streplace(tmp,strfind(":",tmp,n=20),"");
        session_filename = streplace(tmp,strfind(" ",tmp,n=20),"_")+".dat";
      }
      file = session_filename;
    }
    f = open(file,"a");
  }

  if (f) {
    if (coord_name != current_objname) { // write header:
      write,f,format="\n%s\n","============================================";
      write,f,format="%s: %s  %s\n", coord_name, coordX, coordY;
      write,f,format="%s\n","============================================";
      current_objname = coord_name;
    }
    if (n==0) {
      write,f,format="\n%s\n","No valid group found";
      return;
    }
    write,f,format="\nGroup %d\nStrehl avg=%.1f  rms=%.1f  min=%.1f  max=%.1f\n", \
      n,sall(n).avgstrehl*100, sall(n).rmsstrehl*100,                   \
      sall(n).minstrehl*100, sall(n).maxstrehl*100;
  } else {
    write,format="\nGroup %d  (Strehl: avg=%.1f  rms=%.1f  min=%.1f  max=%.1f)\n", \
      n,sall(n).avgstrehl*100, sall(n).rmsstrehl*100,                   \
      sall(n).minstrehl*100, sall(n).maxstrehl*100;
  }

  nstars = numberof(*sall(n).starmag);
  dis = sqrt( (*sall(n).starx)^2. + (*sall(n).stary)^2. );  
  // we want the brightest on CWFS 3:
  wfs_order = [3,2,1];
  
  for (i=1;i<=nstars;i++) {
    // RA & DEC
    ra  = ra_decimal2hms((*sall(n).starra)(i));
    dec = dec_decimal2dms((*sall(n).stardec)(i));
    fmt = "GS%d-%d  %02.0f:%02.0f:%06.3f   %+02.0f:%02.0f:%05.2f   R=%.2f  dist=%.1f\"\n";
    if (f) {
      write,f,format=fmt, wfs_order(i),n,ra(1),ra(2),ra(3),     \
        dec(1),dec(2),dec(3),(*sall(n).starmag)(i),dis(i);
    } else {
      write,format=fmt,wfs_order(i),n,ra(1),ra(2),ra(3),        \
        dec(1),dec(2),dec(3),(*sall(n).starmag)(i),dis(i);
    }
  }
  if (f) {
    write,format="\nAsterism %d saved in %s\n",n,file;
    close,f;
  }
}

