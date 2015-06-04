require,"mascot_utils.i";

func disp_dss_image(void)
{
  window,1;
  fma;
  psize = 0.48;
  sx = dimsof(image_survey_vis)(2)*psize;
  sy = dimsof(image_survey_vis)(3)*psize;

  pli,image_survey_vis,-sx/2.,-sy/2.,sx/2,sy/2.;

  status = overlay_apertures();
}

func disp_2mass_image(void)
{
  window,2;
  fma;
  psize = 0.48;
  sx = dimsof(image_2mass)(2)*psize;
  sy = dimsof(image_2mass)(3)*psize;

  pli,image_2mass,-sx/2.,-sy/2.,sx/2,sy/2.;
  
  status = overlay_apertures();
}


func overlay_apertures(void)
{
  t = span(0.,2*pi,200);
  rad = 60.;
  x = rad*cos(t);
  y = rad*sin(t);

  plg,y,x,color="white";
  hf = halffield;
  plg,[-hf,-hf,hf,hf,-hf],[-hf,hf,hf,-hf,-hf],color="red";
}


func disp_stars(void)
/* DOCUMENT disp_stars(void)
   Display the valid star positions in DSS and 2MASS images
   SEE ALSO:
 */
{
  if (starlist==[]) {
    status = disp_dss_image();
    status = disp_2mass_image();
    return;
  }
  
  ns = dimsof(starlist)(0);

  
  // display vis image (refresh)
  status = disp_dss_image();
  
  // display 2mass image (refresh)
  if (get_2mass_image) {
    status = disp_2mass_image();
    // window,2;
    // fma;
    // pli,image_2mass;
    nbdisp=2;
  } else nbdisp=1;

  for (n=1;n<=nbdisp;n++) {
    window,n;
    for (i=1;i<=ns;i++) {
      // dd = abs(starlist(1,i),starlist(2,i)); // distance to center;
      // if (dd<(1.4*60)) {
      plp,starlist(2,i),starlist(1,i),symbol=12,width=5,size=0.8,color="red";
      // }
    }
  }
}



func disp_prev_next_asterism(dir)
{
  extern current_ast;

  nast = numberof(sall);

  n = current_ast+dir;

  if (n>nast) exit,"Asterism# > available number";
  if (n<1)    exit,"Asterism# <1";

  status = print_asterism(n);

  window,3;
  disp_strehl_map,sall(n);
  
  current_ast = n;

  // display asterism outline in images:
  status = disp_stars();
  x = (*sall(n).starx);  x = _(x,x(1));
  y = (*sall(n).stary);  y = _(y,y(1));
  window,1; plg,y,x,color=char([200,200,200]);
  if (get_2mass_image) { window,2; plg,y,x,color=char([200,200,200]);}

  // set sensitivity of next/previous asterism in GUI:
  if (current_ast==1)    pyk,"glade.get_widget('prev_ast').set_sensitive(0)";
  if (current_ast==nast) pyk,"glade.get_widget('next_ast').set_sensitive(0)";
  if (current_ast>1)     pyk,"glade.get_widget('prev_ast').set_sensitive(1)";
  if (current_ast<nast)  pyk,"glade.get_widget('next_ast').set_sensitive(1)";
}



func disp_strehl_map(sdata)
/* DOCUMENT disp_strehl_map(sdata)
   Display a nice contour plot of the Strehl map, with overlayed
   GeMS entrance aperture, GSAOI FoV, TTGS location.
   SEE ALSO:
 */
{
  local nfp;

  fma;
  
  xloc     = *sdata.starx;
  yloc     = *sdata.stary;
  mag      = *sdata.starmag;
  hf       = sdata.halffield;
  maxtheta = 60.;
  
  dimx     = dimsof(*sdata.strehl_map)(2);
  nfp      = dimx/2;
  tx       = (float(indgen(2*nfp+1))-nfp-1)/nfp*maxtheta;
  ty       = (float(indgen(2*nfp+1))-nfp-1)/nfp*maxtheta;

  nplot = numberof(tx) - 1;
  plfc,transpose(double(*sdata.strehl_map)),ty(,-::nplot),tx(-::nplot,), \
    levs=(indgen(21)-1)/20.;

  for (i=1;i<=21;i++) {
    Text = swrite((i-1)/20.,format="%.2f");
    st = long(((i-1)/20.)*100.);
    plc,transpose(double(*sdata.strehl_map)),ty(,-::nplot),tx(-::nplot,),
      levs=[(i-1)/20,i/20.];
  }

  xytitles,"X Offset [arcsec]","Y Offset [arcsec]",[0.023,0.016];

  plg,_(yloc,yloc(1)),_(xloc,xloc(1)),color=[30,30,30];
  
  plmk,yloc,xloc,marker=6,msize=.3,width=10,color="black";
  plg,[-hf,-hf,hf,hf,-hf],[-hf,hf,hf,-hf,-hf],color="red";
  plg,60.*cos(2*pi*(indgen(100)-1)/99.),60.*sin(2*pi*(indgen(100)-1)/99.),\
    color="blue";

  ns = numberof(mag);
  for (i=1;i<=ns;i++) {
    text = swrite(float(mag(i)),format="%.1f");
    plt,text,xloc(i),yloc(i),justify="CB",tosys=1,color="black";
  }

  maxtheta = 60;
  limits,-maxtheta,maxtheta,-maxtheta,maxtheta;
  colorbar,adjust=-0.017,levs=19;

  // update GUI strehl fields:
  pyk,swrite(format="glade.get_widget('avgstrehl').set_text('%.1f')",sdata.avgstrehl*100);
  pyk,swrite(format="glade.get_widget('rmsstrehl').set_text('%.1f')",sdata.rmsstrehl*100);
  pyk,swrite(format="glade.get_widget('minstrehl').set_text('%.1f')",sdata.minstrehl*100);
  pyk,swrite(format="glade.get_widget('maxstrehl').set_text('%.1f')",sdata.maxstrehl*100);

}

