// mascot Strehl compute/optimize using distortion modes
// instead of quadratic/tt phase (original method).

// todo:
// - introduce wind velocity in zernike spectra (done)
// - make sure ftcb is not sqrt (done, ok)
// - make sure it works with 2 stars and 1 star only (done)
// - do I need to isolate TT from quads? (done: not needed)
// - make the thing faster (compute what takes time)
// - transform Cn2 profile into distortion amplitude (done for
//   tt, approximated for quads).
// - clean up code implementation

// checks
// * when putting wind=10 for all layers, knee for tt is 0.63Hz
// and for quad at 1.10Hz. That matches well with the formula
// 0.3*(n+1)(v/D), which would give 0.3*2*10/8.= 0.75 and 1.12Hz.
// so there seems to be no issue with null_modes_spectra (v>=0.5)
// * changing the number of point npt in over which ftcb is
// computed doesn't change significantly the results (lowering
// from 512 to 256 seems to reduce strehl from 1% to 2 % when @ 60%
// settling on 256
// * changing optim_npt, the number of point on which the strehl is estimated
// for the minimization, does not change the result (tried 5, 11, 21,
// settling on 5).

// Authors, Francois Rigaut,  2010, for this file.
//          Damien Gratadour, 2008-2010.

// many parameters we use here are defined in mascot.conf:
require,"mascot.conf";
require,"mascot_utils.i";
require,"mascot_disp.i";
require,"amoeba.i";

nmodes  = 5;


struct strehl_str
{
  double avgstrehl;
  double rmsstrehl;
  double minstrehl;
  double maxstrehl;
  double halffield;
  pointer strehl_map;
  pointer strehl_map_halffield;
  pointer tiperr;
  pointer tilterr;
  pointer starx;
  pointer stary;
  pointer starra;
  pointer stardec;
  pointer starmag;
}


func create_distortion(mn,x,y)
/* DOCUMENT create_distortion(mn,x,y)
   Returns local distortion at (x,y) for mode N in arcsec
     
   mn = mode n (out of 5)
   x and y = coordinates in the fov (say in arcsec). can be arrays.
   Modes normalized so that modes 1 & 2 (TT) return always 1
   Modes 3,4,5 return 1 for a distance of 100 arcsec for mode3.
   4 and 5 normalized as 3 (that is considering quadratic have same
   variance in turb and taking into account Z4,5,6 norm factor).
   SEE ALSO:
 */
{
  if (mn==1) return [x*0+1.,x*0.];
  else if (mn==2) return [x*0.,x*0+1.];
  else if (mn==3) return [x,y]/100.;
  else if (mn==4) return 2*sqrt(6)/4/sqrt(3.)*[y,x]/100.;
  else if (mn==5) return 2*sqrt(6)/4/sqrt(3.)*[x,-y]/100.;
  else error,"mn out of range";
}


func get_distortion_vfield(mnv,npt,halffield,&xloc,&yloc,offset=,x=,y=,plot=)
/* DOCUMENT get_distortion_vfield(mnv,npt,halffield,&xloc,&yloc,offset=,x=,y=,plot=)
   Return a distortion vector field (from modes with input coefficients) at
   given points in the field of view.
   
   mnv = vector of nmodes mode coefficients, e.g. [1,2,0,0,0]
   x & y = specific coordinates at which to compute the distortions
   if x is not set, then this routine returns the distortion vector field
   at points defined on a square grid, with npt x npt points covering
   a field of view = 2 * halffield.
   if square grid is requested, xloc & yloc on ouput are the grid location
   at which the distortion were computed.
   plot= as it says. Plots the resulting vector field.
   SEE ALSO:
 */
{
  if (x==[]) square=1; // no specific coordinates are supplied (x & y)

  if (offset==[]) offset=[0.,0.];
  
  if (square) { // compute grid locations, covering [-halffield,halffield]
    xy = indices(npt)-(npt+1)/2.;
    xy = xy/max(xy)*halffield-offset(-,-,);
    x = xloc = xy(,,1);
    y = yloc = xy(,,2);
  }

  // get distortion vector field for first mode.
  d = mnv(1)*create_distortion(1,x,y);

  // add other modes.
  for (i=2;i<=numberof(mnv);i++) d += mnv(i)*create_distortion(i,x,y);

  // possible plot if requested.
  if (square && plot) {
    fma; 
    plvf,d(*,2),d(*,1),xy(*,2),xy(*,1),autoscale=1;
    plmargin;
  }
  
  return d;
}


func mascot_compute_strehl(void)
/* DOCUMENT mascot_compute_strehl(void)
   Main routine. originally from Damien Gratadour.
     
   SEE ALSO:
 */
{
  ns = dimsof(starlist)(0);

  gso = array(double,[2,2,ns]);
  mag = ra = dec = array(double,ns);
  for (i=1;i<=ns;i++) {
    gso(,i) = starlist(1:2,i);
    mag(i)  = starlist(5,i);
    ra(i)   = starlist(10,i);
    dec(i)  = starlist(11,i);
  }

  sdata = mascot_optimize(gso,mag);
  sdata.starra  = &ra;
  sdata.stardec = &dec;
  
  return sdata;
}



func mascot_optimize(gso,mag)
/* DOCUMENT mascot_optimize(gso,mag)
   
   Optimizes the modes & mode gains.
   Computes the modal cmat, etc
   Call amoeba that maximizes a criteria (strehl avg/rms)
   Compute final performance: Strehl map
   gso = array(coords,nstars)
   mag = array(nstars)
   SEE ALSO:
*/
{
  extern nmodes_cont, mcmat, mprop;
  // extern ev, mta;
  extern dfields;
  
  nstars = numberof(mag);

  sdata = strehl_str();
  sdata.halffield = halffield;
  sdata.starx     = &gso(1,);
  sdata.stary     = &gso(2,);
  sdata.starmag   = &mag;
  
  // first we recenter the asterism.
  // We do that because we want to keep the segregation between TT and
  // quadratic modes. If we didn't do that, in case we have for instance,
  // 3 stars almost aligned off-axis, then we will necessarily have
  // some TT in the mode with the lowest eigenvalue, which means
  // we will also have some quadratic in the first 2 modes. This has several
  // drawbacks: (1) the turbulence covariance matrix of these eigenmodes is
  // not diagonal anymore. So we have to take into account extra-diagonal terms.
  // (2) each mode being a mix between TT and quadratics, they are not ranked
  // from highest to lowest variance, thus the modal gain optimization will be
  // less efficient (I feel), and (3) it is more difficult to introduce
  // vibrations + windshake.
  // So, we recenter here, do our stuff, and then estimate the perf on a shifted
  // area.
  gso_off = gso(,avg);
  gso = gso - gso_off(,-);
  
  // compute imat:
  imat = array(0.,[2,6,nmodes]);
  for (mn=1;mn<=nmodes;mn++) {
    signal = [];
    for (i=1;i<=nstars;i++) grow,signal,create_distortion(mn,gso(1,i),gso(2,i));
    imat(1:2*nstars,mn) = signal;
  }


  // prepare imat inversion
  ev = SVdec(imat,u,vt);
  nmodes_cont = sum(ev!=0); // number of controlled modes
  mta = transpose(vt);      // modes to act transfer matrix
  if (debug>2) {
    write,format="%s","eigenvalues: "; write,format="%.3g  ",ev; write,"";
  }
  
  // Build modal cmat:
  eem1 = unit(nmodes)*0.;
  for (i=1;i<=nmodes;i++) if (ev(i)>0) eem1(i,i) = 1./ev(i);
  mcmat = eem1(,+) * u(,+);


  // noise propagation
  nca = unit(6)*0.;
  ttnv = wfs_noise(mag,verbose=(debug!=0)); // noise per WFS in arcsec rms
  for (i=1;i<=nstars*2;i++) nca(i,i) = ttnv(i);
  mprop = u(+,) * nca(,+);
  mprop = eem1(,+) * mprop(+,);
  // noise propagation matrix on eigenmodes:
  mprop = mprop(,+) * mprop(,+);
  // it is indeed covariance (diagonal = variance)


  // eigenmodes vector-field display, if needed:
  if (debug>2) {
    window,5;
    tv,mta;
    window,6;
    tv,mprop;
    window,4;
    for (i=1;i<=nmodes_cont;i++) {
      get_distortion_vfield, mta(,i), 10, 40, plot=1;
      plp, gso(2,), gso(1,), symbol=20, color="red";
      pltitle,swrite(format="eigenmode #%d",i);
      hitReturn;
    }
  }
  
  // pre-compute distortion fields:
  dfields = array(0.,[4,nmodes,optim_npt,optim_npt,2]);
  // first index: mode #
  // second and third: spatial X and Y
  // fourth index: local distortion value (1:x or 2:y)
  // so that dfields(2,,,1)
  // is the X component of the distortion field created by 2nd mode.
  // this is a 2D array covering a field = 2*halffield with
  // optim_npt x optim_npt points.
  for (i=1;i<=nmodes;i++) 
    dfields(i,,,) = get_distortion_vfield(mta(,i),optim_npt,halffield,xloc,yloc);

  // call minimization routine. This should return the best gains
  bg = amoeba(0.01, get_strehl_map,nc, fval, nMax=1000,
              p0=-0.7*array(1.,nmodes_cont),scale=0.2);
  
  if (debug>1) {
    write,format="%s ","best gains = "; write,format="%.3f  ",10^bg; write,"";
  }
  

  // now compute the final Strehl map given the gains we just found.
  // pre-compute distortion fields:
  dfields = array(0.,[4,nmodes,smap_npt,smap_npt,2]);
  for (i=1;i<=nmodes;i++) 
    dfields(i,,,) = get_distortion_vfield(mta(,i),smap_npt,halffield,\
                                          xloc,yloc,offset=gso_off);

  // get the strehl map:
  get_strehl_map,bg,smap,tiperr,tilterr;
  sdata.avgstrehl = avg(smap);
  sdata.rmsstrehl = (smap)(*)(rms);
  sdata.minstrehl = min(smap);
  sdata.maxstrehl = max(smap);
  sdata.strehl_map_halffield = &smap;
  
  write,format="Strehl over %.1f\": avg=%.1f  rms=%.1f  min=%.1f  max=%.1f\n",
    halffield*2, sdata.avgstrehl*100, sdata.rmsstrehl*100,
    sdata.minstrehl*100, sdata.maxstrehl*100;

  // Finaly compute the final strehl map, but now over the whole 60"x60" FoV
  // pre-compute distortion fields:
  dfields = array(0.,[4,nmodes,smap_npt,smap_npt,2]);
  for (i=1;i<=nmodes;i++) 
    dfields(i,,,) = get_distortion_vfield(mta(,i),smap_npt,60.,\
                                          xloc,yloc,offset=gso_off);
  
  // get the strehl map:
  get_strehl_map,bg,smap,tiperr,tilterr,stop=(debug>10);
  sdata.strehl_map = &smap;
  sdata.tiperr     = &tiperr;
  sdata.tilterr    = &tilterr;

  return sdata;
}

func get_strehl_map(lgains,&strehl,&tiperr,&tilterr,stop=)
/* DOCUMENT get_strehl_map(lgains,&strehlmap,&tiperr,&tilterr,stop=)
   Given a noise propagation matrix mprop (passed in extern) and
   the modes distortion vector fields dfields (passed in extern),
   compute a strehl map for the mode gains lgains (input). Take
   into account the turbulence residuals and the propagated noise.
   This function is normally called by amoeba (find best gains
   by minimizing a criteria), but for convenience, can also be
   called directly, in which case it returns the strehl maps
   and the tip/tilt errors.
   SEE ALSO:
 */
{
  // extern nois_cov;
  local npt;
  
  turb = nois = g = array(0.,nmodes);
  g(1:nmodes_cont) = 10.^lgains;
  freq  = sp(,1);
  
  // limits upper freq range for spline:
  if (max(spv(,1))>sampfreq) {
    w = where(spv(,1)<sampfreq)(0);
    spv = spv(1:w,);
  }
  freqv = spv(,1);

  rmsvib = array(0.,2);
  rmsvib(1) = sum((*tipvibrms)^2.);
  rmsvib(2) = sum((*tiltvibrms)^2.);
  
  // compute transfer functions for said gains.
  for (i=1;i<=nmodes;i++) {
    // npt = 256;
    npt = 512;
    hs      = array(0.,[2,npt+1,4]);
    hs(2:,) = ftcb(1./sampfreq,0.3e-3,2e-3,g(i),npt);
    // add missing values for freq=0
    if (g(i)==0) hs(1,) = [0.,0.,1.,0.];
    else hs(1,) = [0.,0.,0.,1.];
    // above these are to be applied on PSD.
    // computed in ftcb as h*conj(h), so OK.
    // servolag for turb + vibrations
    herror = clip(spline(hs(,3),hs(,1),freq),0,);
    turb(i) = sum( rmsmodes(i)^2 * sp(,i+1) * herror );
    if (i<=2) {
      // vibrations + windshake:
      // NO ! we can't be sure 2 first modes are Tip and Tilt.
      // let's add this separately. TO BE DONE.
      // YES. now with re-centered asterism, we have clean
      // TT isolation as modes 1 and 2. DONE.
      herror = clip(spline(hs(,3),hs(,1),freqv),0,);
      if (!novibs) {
        turb(i) += sum( rmsvib(i)^2. * spv(,i+1) * herror );
      }
    }
    // noise
    hnoise = hs(,4)/npt;
    nois(i) = sum( hnoise );
  }

  nois = nois * (g>0.); // just to make sure nois=0 if gain=0
  
  // servolag (variance):
  // turb_var = turb * rmsmodes^2.; // compensated rms / mode
  turb_var = turb; // compensated rms / mode

  // noise (variance)
  nois_cov = (nois(-,))*mprop*(nois(,-));
  if (debug>3) {
    window,6;
    tv,nois_cov;
  }
  
  // map of tt error in field of view, turbulence residuals:
  tiperr  = turb_var(+) * (dfields^2.)(+,,,1);
  tilterr = turb_var(+) * (dfields^2.)(+,,,2);

  // noise contribution, covariance:
  for (i=1;i<=nmodes_cont;i++) {
    for (j=1;j<=nmodes_cont;j++) {
      // if (j!=i) continue; << clearly wrong !!!
      tiperr  += nois_cov(i,j) * dfields(i,,,1) * dfields(j,,,1);
      tilterr += nois_cov(i,j) * dfields(i,,,2) * dfields(j,,,2);
    }
  }
  
  tiperr  = sqrt(tiperr); // in arcsec
  tilterr = sqrt(tilterr);

  tiperr_rd = tiperr*4.848e-6*tel_diam*2*pi/(lambdaim*1e-6)/4.;
  tilterr_rd = tilterr*4.848e-6*tel_diam*2*pi/(lambdaim*1e-6)/4.;

  // see strehl_vs_ttrms() below, this is it, fairly good approximation:
  strehl = sqrt(1./(1.+2.*tiperr_rd^2.))*sqrt(1./(1.+2.*tilterr_rd^2.));

  // other previous attempts
  // strehl = sqrt(1./(1.+2.*tiperr_rd))*sqrt(1./(1.+2.*tilterr_rd));
  // strehl = exp(-tiperr_rd^2.-tilterr_rd^2.);  
  // strehl = sqrt(1./(1.+2.*sqrt(tiperr_rd^2.+tilterr_rd^2.)));

  if (stop) error;

  // return -log(avg(strehl)/(0.05+strehl(*)(rms)));
  // return -log(max(strehl)+0.2*avg(strehl));
  return -log(avg(strehl));
}


//================================================================
// test & check functions
//================================================================

func strehl_vs_ttrms(void)
/* DOCUMENT strehl_vs_ttrms(void)
   This is a check routine, to check the formula Strehl = f (tt rms).
   The conclusion is that the correct formulation is:
   strehl_att = sqrt(1./(1.+2*ttrmsx_rd^2.))*sqrt(1./(1.+2*ttrmsy_rd^2.));
   where ttrmsx_rd is the tip variance, as in phase variance = a_2^2
   SEE ALSO:
 */
{
  dim = 256;
  npt = 30;
  
  bdim = 8*dim;
  pup = dist(bdim)<(bdim/2./16.);
  psf = roll(abs(fft(pup,1))^2.,[dim/2,dim/2])(1:dim,1:dim);
  psftf = fft(roll(psf),-1);
  xy = indices(dim)-(dim+1)/2.
  // tv,psf; // fwhm = lambda/D should be 16 pixels.
  satt = array(0.,[2,npt,npt]);
  satt_theo = array(0.,[2,npt,npt]);
  for (i=1;i<=npt;i++) {
    ttrmsx = i*2./npt; // rms in lambda/D units, from 0 to 2.
    fwhmx = ttrmsx * 16 * 2.35;
    for (j=1;j<=npt;j++) {
      ttrmsy = j*2./npt; // rms in lambda/D units, from 0 to 2.
      fwhmy = ttrmsy * 16 * 2.35;
      // hence FWHM of corresponding gaussian from 0 to 2*2.35*16
      g = exp(-xy(,,1)^2./(fwhmx/1.66)^2. - xy(,,2)^2./(fwhmy/1.66)^2.);
      g = g/sum(g);
      im = abs( fft( fft(g,1) * psftf, -1) );
      tv,im;
      satt(i,j) = max(im);

      /* normalization ttrms from lambda/D units to rd of phase:
         ttrmsx_rd  = ttrmsx; // in lambda/D units
         lambda = d = 1.; // will normalize out.
         ttrmsx_rd *= lambda/d; // in rd of angle
         ttrmsx_rd *= d; // in phase difference [m] at edges
         ttrmsx_rd *= (2*pi/lambda); // in phase diff [rd] at edges
         ttrmsx_rd *= (1./4.); // in Z2 coefficients
         that is, in consolidating:
         ttrmsx_rd = ttrmsx * pi/2.
      */
      ttrmsx_rd = ttrmsx * pi/2.;
      ttrmsy_rd = ttrmsy * pi/2.;
      satt_theo(i,j) = sqrt(1./(1.+2*ttrmsx_rd^2.))*sqrt(1./(1.+2*ttrmsy_rd^2.));
    }
  }
  satt = satt/max(satt);

  // ok, this works good. we have satt = satt_theo.
  // so this proves that the approximate theoretical expression
  // of Strehl attenuation vs tt rms is indeed:
  // strehl_att = sqrt(1./(1.+2*ttrmsx_rd^2.))*sqrt(1./(1.+2*ttrmsy_rd^2.));
  // where ttrmsx_rd is the tip variance, as in phase variance = a_2^2
  
  return [satt,satt_theo];
}

func snr_vs_nph_quadcell(&nph_star,&nph_sky,nit=)
/* DOCUMENT snr_vs_nph_quadcell(void)
   Check theoretical formulation of noise against simulation,
   for a quadcell.
   Note the expression of denominator in quadcell formula. I divide
   by useful # of photons, excluding dark and sky.
   
   SEE ALSO:
 */
{
  // total # ph/frame
  nph_star = 10^((indgen(7)-1.)/2.)(::-1);
  // sky ph#/channel/frame, including dark:
  nph_sky  = _(0.,10^((indgen(5)-3.)/2.));
  write,format="%s","nph_star = "; write,format="%.2f  ",nph_star; write,"";
  write,format="%s","nph_sky  = "; write,format="%.2f  ",nph_sky; write,"";
  if (!nit) nit = 100000;
  res      = array(0.,[3,numberof(nph_star),numberof(nph_sky),2]);
  sig      = array(0.,nit);
  sigok    = array(0,nit);
  // 1rst plan = rms
  // 2nd plan  = signal.
  write,format="%s","nph_star = ";
  for (nst=1;nst<=numberof(nph_star);nst++) {
    write,format="%.2f ",nph_star(nst);
    sig *= 0.; sigok *= 0;
    for (nsk=1;nsk<=numberof(nph_sky);nsk++) {
      for (i=1;i<=nit;i++) {
        phst = poidev( array(nph_star(nst)/4.,4) );
        // phst = array(0.,4);
        // phst(1:2) = poidev( array(nph_star(nst)/2.*0.6) );
        // phst(3:4) = poidev( array(nph_star(nst)/2.*0.4) );
        phsk = poidev( array(nph_sky(nsk),4) );
        ph = phst + phsk;
        sigok(i) = (sum(phst)>0);
        if (sum(phst)>0) sig(i)=(ph(1:2)(sum)-ph(3:4)(sum))/sum(phst);
      }
      w = where(sigok);
      res(nst,nsk,1) = sig(w)(rms);
      res(nst,nsk,2) = avg(sig(w));
    }
  }
  // write,format="\nminmax(signal) (should be 0.2) = %.3f  %.3f\n",min(res(,,2)),max(res(,,2));
  write,format="\nminmax(signal) (should be 0) = %.3f  %.3f\n",min(res(,,2)),max(res(,,2));
  
  return res;
}
