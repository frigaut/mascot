/*********************************************/
/* Python-Yorick interface related functions */
/*********************************************/

func pyk_status_push(msg,id=)
{
  if (id==[]) id=1;
  pyk,swrite(format="pyk_status_push(%d,'%s')",id,msg);
}

func pyk_status_pop(id=)
{
  if (id==[]) id=1;
  pyk,swrite(format="pyk_status_pop(%d)",id);
}

func pyk_info(msg)
{
  if (numberof(msg)>1) msg=sum(msg+"\\n");
  // or streplace(msg,strfind("\n",msg),"\\n")
  pyk,swrite(format="pyk_info('%s')",msg)
}

func pyk_info_w_markup(msg)
{
  if (numberof(msg)>1) msg=sum(msg+"\\n");
  // or streplace(msg,strfind("\n",msg),"\\n")
  pyk,swrite(format="pyk_info_w_markup('%s')",msg);
}

func pyk_error(msg)
{
  if (numberof(msg)>1) msg=sum(msg+"\\n");
  // ok, here the problem is that "fatal errors", when called from shell,
  // should bail you out (quit yorick). But if they do, then the python
  // process is also killed and then the error message never appears on screen.
  // thus the use of zenity in *all* cases.
  //  if (_pyk_proc) {
  //    pyk,swrite(format="pyk_error('%s')",msg);
  //  } else { // python not started yet, use zenity
    system,swrite(format="zenity --error --text=\"%s\"",msg);
    //  }
}

func pyk_warning(msg)
{
  if (numberof(msg)>1) msg=sum(msg+"\\n");
  pyk,swrite(format="pyk_warning('%s')",msg)
}

func gui_progressbar_frac(frac) {
  frac = clip(frac,0.,1.);
  pyk,swrite(format="progressbar.set_fraction(%f)",float(frac));
}


func gui_progressbar_text(text) {
  pyk,swrite(format="progressbar.set_text('%s')",text);
}


func gui_message(msg)
{
  pyk,swrite(format="statusbar.push(1,'%s')",msg);
}


/*********************************************/
/* Python-Yorick interface related functions */
/*********************************************/

func mascot_win_init(pid1,pid2,pid3,pid4)  // main display window
{
  //  xft,0;
  wp = (window_pad!=[]?window_pad:0);
  // DSS image:
  // window,1,dpi=66,wait=1,parent=pid1,ypos=-25,style="nobox.gs";
  window,1,dpi=63,wait=1,parent=pid1,ypos=-3+wp,xpos=-3+wp,style="nobox.gs";
  // require,"png.i";
  // logo = png_read("gems_logo_nobkg.png");
  // pli,logo(1:3,,::-1);
  limits,square=1;
  // pli,dist(300);
  // limits,square=1;

  // 2MASS image
  window,2,dpi=63,wait=1,parent=pid2,ypos=-3+wp,xpos=-3+wp,style="nobox.gs";
  limits,square=1;

  // Strehl map
  window,3,dpi=63,wait=1,parent=pid3,ypos=-25+wp,xpos=-3+wp,style="pymex.gs";

  // zoom
  //  window,4,dpi=41,wait=1,parent=pid4,ypos=-25,style="pymex.gs";

  window,1;
  pyk,"done_init = 1";

  pyk,swrite(format="glade.get_widget('avg_rms_crit').set_value(%f)",
             float(avg_rms_criteria));
  pyk,swrite(format="glade.get_widget('minmag').set_text('%.1f')",
             float(mag_min_threshold));
  pyk,swrite(format="glade.get_widget('maxmag').set_text('%.1f')",
             float(mag_max_threshold));
  pyk,swrite(format="glade.get_widget('nstar_limit').set_text('%d')",
             long(nstar_limit));
  // on myst I used this with a simple .gtkrc to force using clearlook
  // in remote display:"
  // pyk,"glade.get_widget('skin2').activate()";

  msg = swrite(format="Mascot %s ready",mascot_version);
  write,msg;
  maybe_prompt;
  pyk_status_push,msg;
  
  if (obj2process!=[]) process_multi,obj2process;
  
  //  mascot_lut,0;
}

func delete_star_under_cursor(void)
{
  extern starlist;
  
  if (starlist==[]) {
    msg = "Star list empty !";
    write,msg;
    pyk_status_push,msg;
    return; // no stars available
  }
  
  cw = current_window();
  if (cw<0) return; // there's no window

  m = current_mouse();
  if (m==[]) return; // not currently over a window.
  
  m = m(1:2);
  d = sqrt( (starlist(1,)-m(1))^2. + (starlist(2,)-m(2))^2. );
  w = where(d!=min(d));
  starlist = starlist(,w);
  status = disp_stars();
}

func info_star_under_cursor(void)
{
  if (starlist==[]) {
    msg = "Star list empty !";
    write,msg;
    pyk_status_push,msg;
    return; // no stars available
  }
  
  cw = current_window();
  if (cw<0) return; // there's no window

  m = current_mouse();
  if (m==[]) return; // not currently over a window.
  
  m = m(1:2);
  d = sqrt( (starlist(1,)-m(1))^2. + (starlist(2,)-m(2))^2. );
  w = where(d==min(d))(1);
  msg = swrite(format="Star offsets = (%.1f,%.1f) R=%.2f",
               starlist(1,w),starlist(2,w),starlist(5,w));
  write,msg;
  pyk_status_push,msg;
}


func recenter_here(void)
{
  cw = current_window();
  m = current_mouse()(1:2); // offsets to perform

  dec = decdec = dec_dms2decimal(ytarget);
  dec += m(2)/3600.;
  dec = dec_decimal2dms(dec);

  ra = ra_hms2decimal(xtarget);
  ra -= m(1)/3600./cos(decdec*pi/180.);
  ra = ra_decimal2hms(ra);
  
  cx = swrite(format="%02.0f:%02.0f:%06.3f",ra(1),ra(2),ra(3));
  cy = swrite(format="%+02.0f:%02.0f:%05.2f",dec(1),dec(2),dec(3));

  process,cx+" "+cy;
}


func process_multi(objects)
{
  extern object_list;
  extern autoflag;
  
  if (objects!=[]) {
    object_list = objects; // init.
    // autoflag = 1;
  } else {
    // this is at least the second object.
    // we're going to save the first object now, as it is
    // a convenient, if not logical, place to do so.
    nast = numberof(sall);
    if ( (autoflag) && (nast==0) ) print_asterism,0,1;
    if ( (autoflag) && (nast>0) ) {
      for (i=1;i<=min([nast,nast2print]);i++) print_asterism,i,1;
    }
  }

  if (object_list == []) return; // nothing left to process. we're done

  // put next object in queue for processing:
  object=object_list(1);

  // get rid of first object in list:
  if (numberof(object_list)>1) object_list = object_list(2:);
  else object_list=[];

  if (mascot_user_function) status = mascot_user_function(first=(objects!=[]));
  
  process,object;
}


func process(object)
{
  pyk,swrite(format="glade.get_widget('targetEntry').set_text('%s')",object);
  pyk,"glade.get_widget('targetEntry').activate()";
}

func mascot_lut(lut)
{
  require,"idl-colors.i";
  extern rlut,glut,blut;
  local r,g,b;
  extern active_window;

  if (lut!=[]) {  // then read and set new lut
    if (lut==0) palette,"earth.gp";
    else loadct,lut;
    palette,query=1,rlut,glut,blut;  // store
  }

  // invert?
  if (mascot_invertlut) {
    r=rlut(::-1); g=glut(::-1); b=blut(::-1);
  } else {
    r=rlut; g=glut; b=blut;
  }

  // itt:
  if (mascot_itt==1) { // linear
    ind = span(0.,1.,mascot_ncolors);
  } else if (mascot_itt==2) { // sqrt
    ind = sqrt(span(0.,1.,mascot_ncolors));
  } else if (mascot_itt==3) { // square
    ind = (span(0.,1.,mascot_ncolors))^2.;
  } else if (mascot_itt==4) { // log
    ind = log10(span(10.^(-mascot_log_itt_dex),1.,mascot_ncolors)); // 8 dex
    ind -= min(ind);
    ind /= max(ind);
  }
  ind = round(ind*(mascot_ncolors-1)+1);
  r = r(ind); g = g(ind); b = b(ind);

  // and finally, load the palette:
  window,active_window;
  palette,r,g,b;
}


func unzoom_all(void)
{
  for (i=1;i<=3;i++) {
    window,i;
    unzoom;
  }
  window,1;
}

func get_cursor(void)
/* DOCUMENT get_cursor()
   returns [xpos,ypos,system]
   where xpos,ypos is the current pixel over which the cursor are, system 1
   Returns [] if not in correct window;
   SEE ALSO:
 */
{
  extern zoom_window;

  cur = current_mouse(long(zoom_window));
  if (cur==[]) return;
  cur = long(cur)
  cur(1:2) = cur(1:2)+1; //ceil
  return cur;
}

func rad4zoom_incr(void)
{
  rad4zoom=min(rad4zoom+1,mascot_dims(2)/2);
}

func rad4zoom_decr(void)
{
  rad4zoom=max(rad4zoom-1,0);
}

func mascot_zoom
{
  extern from_disp,stop_zoom;
  extern prevxy;
  extern image_survey_vis,image_2mass;
  extern zoom_window;

  if (zoom_window == 1) {
    mascot_im = image_survey_vis;
    window,1;
  }
  if (zoom_window == 2) {
    if (get_2mass_image==0) return;
    mascot_im = image_2mass;
    window,2;
  }
  if (zoom_window == 3) {
    return;
    mascot_im = imageStrehl;
    window,3;
  }

  if (stop_zoom) {
    stop_zoom=0;
    return;
  }

  cur = get_cursor();
  if ( (cur==[]) || (from_disp==4) ) { // not in correct window
    after,0.05,mascot_zoom;
    return;
  }

  if (allof(prevxy==cur(1:2)) && (prevz==rad4zoom)) {  // same positon as before
    after,0.05,mascot_zoom;
    return;
  }

  sys = cur(3);
  if (sys!=0) {

    i = clip(cur(1),1,mascot_dims(2));
    j = clip(cur(2),1,mascot_dims(3));
    //    pyk,swrite(format="y_set_xyz('%d','%d','%4.7g')", \
    //              i,j,float(mascot_im(i,j)));

    local_rad=5;
    x1 = clip(i-local_rad,1,mascot_dims(2));
    x2 = clip(i+local_rad,1,mascot_dims(2));
    y1 = clip(j-local_rad,1,mascot_dims(3));
    y2 = clip(j+local_rad,1,mascot_dims(3));
    //    pyk,swrite(format="y_text_parm_update('localmax','%4.7g')",   \
    //              float(max(mascot_im(x1:x2,y1:y2))));
    sim=mascot_im(x1:x2,y1:y2);
    wm = where2(sim==max(sim))(,1)-local_rad-1;

    // FIXME: blue cursor is not at correct position when pointing lower than smaller indice

    x1 = i-rad4zoom;
    x2 = i+rad4zoom;
    y1 = j-rad4zoom;
    y2 = j+rad4zoom;
    if (x1<1) { x1=1; x2=x1+(2*rad4zoom+1); }
    else if (x2>mascot_dims(2)) { x2=mascot_dims(2); x1=x2-(2*rad4zoom+1); }
    if (y1<1) { y1=1; y2=y1+(2*rad4zoom+1); }
    else if (y2>mascot_dims(3)) { y2=mascot_dims(3); y1=y2-(2*rad4zoom+1); }
    window,4;
    fma;
    if (zoom_cmincmax) pli,bytscl(mascot_im(x1:x2,y1:y2),cmin=cmin,cmax=cmax);
    else pli,mascot_im(x1:x2,y1:y2);
    limits,1,2*rad4zoom,1,2*rad4zoom;
    // plot local maxima location
    //plp,j-y1+wm(2)+0.5,i-x1+wm(1)+0.5,symbol=2,color="blue",width=3;
    // plot cursor location
    plp,j-y1+0.5,i-x1+0.5,symbol=2,color="red",width=3;

    // x and y cuts
    //if (xcut) {
    //  plot_xcut,j;
    //} else if (ycut) {
    //  plot_ycut,i;
    //}

    window,1;
  }
  prevxy = cur(1:2);
  prevz = rad4zoom;
  after,0.05,mascot_zoom;
}

func clear_windows(void)
{
  window,1; fma; redraw;
  window,2; fma; redraw;
  window,3; fma; redraw;
}

func mascot_clean(void)
{
  extern stop_zoom;

  winkill,1;
  winkill,2;
  winkill,3;
  winkill,4;
  stop_zoom=1;
}


func clear_cache(void)
{
  // FIXME: define cache path
  system,"rm cache/*.fits cache/*.dat";
}

func mascot_shortcut_help(void)
{
  help_text = ["The following shortcuts are available:",
               "i : Print info about nearest star",
               "d : Delete nearest star from star list",
               "c or o : Recenter object at mouse (new query)",
               "? or Ctrl+H : This help"];
  write,format="%s\n",help_text;
  pyk_info,help_text;
}
