/* mascot.i
   The MCAO Asterism Configuration Tool.
   syntax: yorick -i mascot.i
   the environment variable MASCOTTOP has to point to the package top directory
   Authors: D. Gratadour, F. Rigaut July 2007

   changelog: see svn logs.
*/

mascot_version = "0.8.2";

require,"pyk.i";
require,"astro_util1.i";
require,"spydr_plugins.i";
require,"pathfun.i";
require,"mascot.conf";
require,"util_fr.i";
require,"histo.i";
require,"plot.i";
require,"mascot_gui.i";
require,"mascot_strehl.i";
require,"mascot_disp.i";

sp = null_modes_spectra();
spv = vib_spectra();

if (get_argv()(0) != "mascot.i") {
  w = where(get_argv()=="mascot.i")(0);
  if (anyof(get_argv()=="@--auto")) {
    w = where(get_argv()=="@--auto")(0);
    autoflag = 1;
  } else autoflag = 0;
  // we'll have to re-parse the args with the introduced
  // separator "@" by the mascot launcher (bash script):
  tmp = get_argv()(w+1:);
  tmp = sum(tmp+" ");
  tmp = strtok(tmp,"@",100);
  tmp = tmp(where(tmp));
  obj2process = strtrim(tmp);
  if (numberof(obj2process)>1) autoflag=1;
 }

/*
  mascot utilities:
  - get target coordinates
  - get star list (positions + magnitudes)
  - get DSS and 2MASS images
 */


func check_mag(bMag,vMag)
/* DOCUMENT check_mag(bMag,vMag)
   Added by FR: I believe this returns the R mag from B and V
   if R is not present.
   SEE ALSO:
 */
{

  a = [0.0530768,0.794274,0.212565,-0.867596,0.77699,-0.161851];
  x = bMag - vMag;

  funcMag = a(1)+a(2)*x+a(3)*x^2.+a(4)*x^3+a(5)*x^4+a(6)*x^5;

  rMag = vMag - funcMag;

  return rMag;
}


func callback1(msg)
/* DOCUMENT callback1(msg)
   Function called back from the shell that takes care of
   fetching the object coordinates.
   The command sent to this shell includes the vizquery for the
   object coordinates and then "over" is echoed.
   When this callback gets the over, it will call the function
   to read and store the object coordinates, and the fetch
   the DSS image from there.
   SEE ALSO:
 */
{
  extern msg1;
  
  grow,msg1,msg;
  if (msg == "over\n") analyse_msg1;
}

func callback2(msg)
/* DOCUMENT callback2(msg)
   Function called back from the shell that takes care of
   fetching the DSS image.
   The command sent to this shell includes the vizquery for the
   DSS image (as fits) and then "over" is echoed.
   When this callback gets the over, it will call the function
   to read, store and display the DSS image.
   SEE ALSO:
 */
{
  extern image_survey_vis;
  extern ytarget;

  if (msg == "over\n") {
    image_survey_vis = fits_read("cache/"+coord_name+"_dss.fits");
    // image_survey_vis = array(0.,[2,300,300]); // FIXME TEMP STATS
    if (numberof(image_survey_vis) != long(300)*long(300)) {
      gui_message,"Could not retrieve the DSS2 image";
    } else {
      status = disp_dss_image();
    }
    if (get_2mass_image) {
      status = get2MASSimage();
    } else {
      window,2;
      plt,"2mass not requested",150,150,tosys=1,justify="CH",height=30;
      pyk,"glade.get_widget('SearchTargetButton').set_sensitive(1)";
      pyk,"glade.get_widget('targetEntry').set_sensitive(1)";
      pyk,"glade.get_widget('SearchStarsButton').set_sensitive(1)";
      pyk,"glade.get_widget('star_filters').set_sensitive(1)";
      if (autoflag) pyk,"glade.get_widget('SearchStarsButton').pressed()";
    }
    gui_message,"Done !";
  }
}

func callback3(msg)
/* DOCUMENT callback3(msg)
   Callback to the call to fetch the star list from a remote
   query to the nomad1 catalog.
   SEE ALSO:
 */
{
  extern msg2;

  grow,msg2,msg;
  if (msg == "over\n") {
    gui_message,"Done !";
    analyse_msg2;
  }
}

func callback4(msg)
/* DOCUMENT callback4(msg)
   Callback to the call to fetch the 2MASS image.
   SEE ALSO:
 */
{
  extern image_2mass;
  extern ytarget;

  if (msg == "over\n") {
    image_2mass = fits_read("cache/"+coord_name+"_2mass.fits")
      if (numberof(image_2mass) != long(300)*long(300)) {
        gui_message,"Could not retrieve 2MASS image";
      } else {
        status = disp_2mass_image();
      }
    window,1;
    pyk,"glade.get_widget('SearchTargetButton').set_sensitive(1)";
    pyk,"glade.get_widget('targetEntry').set_sensitive(1)";
    pyk,"glade.get_widget('SearchStarsButton').set_sensitive(1)";
    pyk,"glade.get_widget('star_filters').set_sensitive(1)";
    if (autoflag) pyk,"glade.get_widget('SearchStarsButton').pressed()";
  }
}

func analyse_msg1(void)
/* DOCUMENT analyse_msg1(void)
   Parse msg1 (answer from remote query to get coordinates for a given
   target name) and store results in coordX/Y (coord string) and
   individual hh/mm/ss etc fields (float)
   SEE ALSO:
 */
{
  extern msg1;
  extern coordX,coordY;
  extern xtarget,ytarget;

  d=strtok(msg1,"\n",100);
  d=d(where(d));
  d=strtok(d(-1),"\t",20);
  coordX=d(5);
  coordY=d(6);
  gui_message,"Target coordinates : "+coordX+" "+coordY;

  a=0.;
  b=0.;
  c=0.;

  sread,coordX,a,b,c;
  xtarget = [float(a),float(b),float(c)];
  sread,coordY,a,b,c;
  ytarget = [float(a),float(b),float(c)];

  gui_message,"Done!";
  write,xtarget;
  write,ytarget;

  if ( allof(xtarget==[0,0,0]) && allof(ytarget==[0,0,0]) ) {
    msg = swrite(format="Can not find target %s in CDS Vizir",coord_name);
    write,msg;
    pyk_status_push,msg;
    pyk,"glade.get_widget('SearchTargetButton').set_sensitive(1)";
    pyk,"glade.get_widget('targetEntry').set_sensitive(1)";
    after,0.0,process_multi;
    return;
  }
                                                  
  
  gui_message,"Looking for DSS image ...";
  
  if (fileExist("cache/"+coord_name+"_dss.fits")) {
    write,format="Image cache/%s_dss.fits found in cache\n",coord_name;
    command="echo over\n";
  } else {
    command = "wget --output-document=cache/"+coord_name+               \
      "_dss.fits 'skyview.gsfc.nasa.gov/cgi-bin/images?Survey=DSS2R&position="+ \
      coordX+", "+coordY+"&size=0.04&Return=FITS'"+"\necho over\n";
  }
  // command = "echo over\n"; // FIXME TEMP STATS
  
  tcsh2,command;

}

func analyse_msg2(void)
/* DOCUMENT analyse_msg2(void)
   Parse the nomad data file (star list) dumped by the findnomad1 call
   and stuff results (GS coordinates and magnitudes) in extern variables
   (mostly starlist).
   SEE ALSO:
 */
{
  extern msg2;
  extern nStars,starlist,nStarsUpdated;
  extern allstarlist;
  extern xtarget,ytarget;

  gui_message,"Analyzing Stars";

  a = b = c = 0.;

  xtargetDeg = (xtarget(1)+xtarget(2)/60.+xtarget(3)/3600.)*15.;
  if (ytarget(1) <0.) ytargetDeg = (ytarget(1)-ytarget(2)/60.-ytarget(3)/3600.);
  else ytargetDeg = (ytarget(1)+ytarget(2)/60.+ytarget(3)/3600.);

  file = rdfile("cache/"+coord_name+"_nomad1.dat");
  write,file;
  d=strtok(file,"\n",1000);
  d=d(where(d));
  fin = where(strmatch(d,"matches"))(0);
  if (fin==6) goto nostar; // nomad file contains no stars
  d=d(6:fin-1);
  nStars = numberof(d);
  nStarsUpdated = nStars;
  // need to be double to maintain accuracy for coordinates (fields 10 & 11)
  starlist=array(double,[2,11,nStars]);
  // field 1 = RA offset
  // field 2 = dec offset
  // field 5 = R mag
  // field 10 = RA
  // field 11 = dec

  for (cptStars=1;cptStars<=nStars;cptStars++) {
    tab1=strtok(d(cptStars),"|",9);
    if (strmatch(tab1(3),"-")) {
      tab2=strtok(tab1(3),"-",2)(1);
      sread,tab2,a;
      xStar = a;
      distStarX = - xStar + xtargetDeg;
      distStarX *= 3600.;
      tab3=strtok(strtok(tab1(3),"-",2)(2)," ",2)(1);
      sread,tab3,a;
      yStar = -1.*a;
      decRad = yStar*pi/180.;
      distStarY = yStar - ytargetDeg;
      distStarY *= 3600.;
      starlist(2,cptStars) = distStarY;
      starlist(1,cptStars) = distStarX*cos(decRad);
    } else {
      tab2=strtok(tab1(3),"+",2)(1);
      sread,tab2,a;
      xStar = a;
      distStarX = - xStar + xtargetDeg;
      distStarX *= 3600.;
      tab3=strtok(strtok(tab1(3),"+",2)(2)," ",2)(1);
      sread,tab3,a;
      yStar = a;
      decRad = yStar*pi/180.;
      distStarY = yStar - ytargetDeg;
      distStarY *= 3600.;
      starlist(2,cptStars) = distStarY;
      starlist(1,cptStars) = distStarX*cos(decRad);
    }
    
    sread,tab2,a;
    starlist(10,cptStars) = a;
    starlist(11,cptStars) = yStar;


    tab2=strtok(tab1(6)," ",4)(1:3);
    for (j = 1;j<=3;j++) {
      if (tab2(j) != "---") {
        sread,tab2(j),a;
        starlist(2+j,cptStars) = a;
      } else starlist(2+j,cptStars) = -27.0;
    }

    tab2=strtok(tab1(7)," ",4)(1:3);
    for (j = 1;j<=3;j++) {
      if (tab2(j) != "---") {
        sread,tab2(j),a;
        starlist(5+j,cptStars) = a;
      } else starlist(5+j,cptStars) = -27.0;
    }

    if (starlist(5,cptStars) == -27) {
      if ((starlist(4,cptStars) != -27) && (starlist(3,cptStars) != -27)) {
        starlist(5,cptStars) = check_mag(starlist(3,cptStars),starlist(4,cptStars));
        if (starlist(5,cptStars) < 18.)
          starlist(9,cptStars) = 1;
        else starlist(9,cptStars) = 0;
      } else {
        if ((starlist(4,cptStars) != -27) && (starlist(4,cptStars) < 17)) starlist(9,cptStars) = 1;
        else starlist(9,cptStars) = 0;
        if ((starlist(3,cptStars) != -27) && (starlist(3,cptStars) < 17)) starlist(9,cptStars) = 1;
        else starlist(9,cptStars) = 0;
      }
    } else {
      if (starlist(5,cptStars) < 18.) starlist(9,cptStars) = 2;
      else starlist(9,cptStars) = 0;
    }

  }

  w = where((abs(starlist(2,)) <= ttgs_max_fov_radius*60.) &
            (abs(starlist(1,)) <= ttgs_max_fov_radius*60.));

  if (numberof(w)==0) { starlist=[]; goto nostar;}

  starlist=starlist(,w);

  // added oct2010 by FR:
  mag = starlist(5,);
  w = where(mag>0);
  if (numberof(w)==0) { starlist=[]; goto nostar;}

  starlist = starlist(,w);

  mag = starlist(5,);
  starlist = starlist(,sort(mag));

  allstarlist = starlist;

  // additionnal magnitude selection:
  status = select_stars_on_mag();

 nostar:
  if (starlist==[]) nstars = 0;
  else {
    nstars = dimsof(starlist)(0);
    pyk,"glade.get_widget('bestAsterismButton').set_sensitive(1)";
    pyk,"glade.get_widget('asterism_criteria').set_sensitive(1)";
  }
  
  pyk,"glade.get_widget('SearchStarsButton').set_sensitive(1)";
  pyk,"glade.get_widget('star_filters').set_sensitive(1)";
  pyk_status_push,swrite(format="Found %d stars with current filter",nstars);

  if (autoflag) pyk,"glade.get_widget('bestAsterismButton').pressed()";
}


func get2MASSimage(void)
/* DOCUMENT get2MASSimage(void)
   Upper command to trigger remote query to get 2MASS image.
   This will trigger spawning the request to a shell, with a
   callback function that will eventually read out the store
   the image when the query is done. The yorick prompt remains
   active during the process.
   SEE ALSO:
 */
{
  extern xtarget,ytarget;

  gui_message,"Looking for 2MASS image ...";
  if (fileExist("cache/"+coord_name+"_2mass.fits")) {
    write,format="Image cache/%s_2mass.fits found in cache\n",coord_name;
    command="echo over\n";
  } else {
    command = "wget --output-document=cache/"+coord_name+\
      "_2mass.fits 'skyview.gsfc.nasa.gov/cgi-bin/images?Survey=2MASSK&position="+\
      coordX+", "+coordY+"&size=0.04&Return=FITS'"+"\necho over\n";
  }
  
  tcsh4,command;

}

func select_stars_on_mag(void)
/* DOCUMENT select_stars_on_mag(void)
   Downselect stars within magnitude range in starlist.
   SEE ALSO:select_stars_not_too_close
 */
{
  extern starlist;
  
  mag = allstarlist(5,);
  w = where( (mag>=mag_min_threshold) & (mag<=mag_max_threshold) );
  starlist = allstarlist(,w);
  
  status = select_stars_not_too_close()
}

func select_stars_not_too_close(void)
/* DOCUMENT select_stars_not_too_close(void)
   Remove from starlist the faint stars around bright stars.
   Here is how it is done:
   starlist is already sorted from the brightest to the faintest star.
   The list is walked, starting from the brightest star.
   For each star, we remove from the starlist all the fainter stars
   that are closer than crowding_radius. Eventually, this provides a list
   of the brightest stars that pad the field as regularly as possible.
   We do that iteratively (increasing crowding_radius at each iteration)
   to end up with no more than nstar_limit stars.
   SEE ALSO:
 */
{
  extern starlist;

  if (starlist==[]) return;
  
  ns = dimsof(starlist)(0);
  valid = array(1,ns);

  crowd_rad = float(crowding_radius);
  
  do {
    for (i=1;i<=ns-1;i++) { // for each stars:
      // look at the distance to next (fainter) stars:
      dd = abs(starlist(1,)-starlist(1,i),starlist(2,)-starlist(2,i));
      ok = (dd>=crowd_rad);
      valid(i+1:) *= ok(i+1:);
    }
    crowd_rad += 2;
  } while (sum(valid)>nstar_limit);
  crowd_rad -= 2;

  write,format="Select stars: found optimum crowding radius=%.0f\"\n",
    crowd_rad;
  
  starlist = starlist(,where(valid));
  
  status = disp_stars();
}

func mascot_search_stars
/* DOCUMENT mascot_search_stars
   Launch the remote query to get the GS list from the nomad catalog.
   SEE ALSO:
 */
{
  window,3;
  fma;

  redraw;
  if (mascot_starcat == 1) {
    gui_message,"Looking for Stars in the Nomad1 Catalog";
    
    if (fileExist("cache/"+coord_name+"_nomad1.dat")) {
      write,format="Image cache/%s_nomad1.dat found in cache\n",coord_name;
      command="echo over\n";
    } else {
      command = cds_bin_path+"findnomad1 -c "+coordX+" "+coordY+       \
        " -r "+swrite(format="%.2f",ttgs_max_fov_radius)+" -E > cache/"+ \
        coord_name+"_nomad1.dat\necho over\n";
    }
    tcsh3,command;
  }
}

func mascot_search_target(object)
/* DOCUMENT mascot_search_target(object)
   Wrapper for the target search function by name (mascot_searchTargetName)
   or by coordinates (mascot_searchTargetCoord)
     
   SEE ALSO:
 */
{
  extern starlist, nast, sall;

  if (strlen(strtrim(object))==0) return;
  
  starlist = [];
  nast = 0;
  sall = [];
  status = clear_windows();

  if (strgrep("[a-zA-Z]",object)(0)!=-1) {
    // this is most likely an object name
    mascot_searchTargetName,object;
  } else {
    // this is most likely object coordinates
    mascot_searchTargetCoord,object;
  }
}

func mascot_searchTargetCoord(coordXY)
{
  extern coordX, coordY;
  extern xtarget, ytarget;
  extern coord_name;

  if (strmatch(coordXY,":")) coordXY=streplace(coordXY,strfind(":",coordXY,n=6)," ");
  tmp = strtok(coordXY," ",6);
  coordX = strtrim(sum(tmp(1:3)+" "));
  coordY = strtrim(sum(tmp(4:6)+" "));

  tmp = strtrim(coordXY);
  coord_name = streplace(tmp,strfind(" ",tmp,n=20),"_");
  
  if ( (strpart(coordY,1:1) != "-") && (strpart(coordY,1:1) != "+") )
    coordY = "+"+coordY;

  a=0.;  b=0.;  c=0.;

  sread,coordX,a,b,c;
  xtarget = [float(a),float(b),float(c)];
  sread,coordY,a,b,c;
  ytarget = [float(a),float(b),float(c)];

  gui_message,"Done!";
  write,xtarget;
  write,ytarget;
  gui_message,"Looking for DSS image ...";
  if (fileExist("cache/"+coord_name+"_dss.fits")) {
    write,format="Image cache/%s_dss.fits found in cache\n",coord_name;
    command="echo over\n";
  } else {
    command = "wget --output-document=cache/"+coord_name+               \
      "_dss.fits 'skyview.gsfc.nasa.gov/cgi-bin/images?Survey=DSS2R&position="+ \
      coordX+", "+coordY+"&size=0.04&Return=FITS'"+"\necho over\n";
  }
  // command = "echo over\n"; // FIXME TEMP STATS
  
  tcsh2,command;

}


func mascot_searchTargetName(name)
{
  extern msg1,msg2,mascot_cat;
  extern coord_name;
  
  coord_name = strtrim(name);
  
  msg1=[];
  msg2=[];

  if (mascot_cat == 1) {
    gui_message,"Looking for coordinates in the History and Accurate Positions for the NGC/IC Objects";
    tcsh1,cds_bin_path+"vizquery -mime=tsv <<====End\n";
    command = "-source=VII/239A/icpos\n-out.form=mini\n-c="+name+"\n====End\necho over\n";
    tcsh1,command;
  }
}

func mascot_select_stars
{
  extern starlistOrig,doneWithSelection,nStarsUpdated;
  extern mascot_pixSize;
  extern image_survey_vis,image_2mass;
  extern selectedIndex;

  npass = 0;
  not_over = 1;
  threshDetect = 2.;
  selectedIndex = [];

  window,1;
  gui_message,"Select stars with mouse, R-click when done, M-click to remove last entry";
  starlistOrig = starlist;

  ndims = dimsof(starlist);
  starlist = array(float,[2,ndims(2),1]);

  do {
    res  = mouse(1,0,"");
    //pyk_status_push,"Processing...";
    c    = res(1:2);
    but  = res(10);
    if (but == 3) not_over = 0;
    if (but == 2) {
      if (npass == 0) {
        gui_message,"You can only unbuffer after having buffered at least one star!";
        pyk_warning,"You can only unbuffer after having buffered at least one star!";
        continue;
      }
      starlist = starlist(,:-1);
      gui_message,"Last measurement taken out of star list";
      //pyk_warning,"Last measurement taken out of star list";
      npass -= 1;
      continue;
    }
    npass += 1;
    if (but == 1) {
      coordX = (c(1)-149)*mascot_pixSize;
      coordY = (c(2)-149)*mascot_pixSize;
      gui_message,swrite(format="%f %f",coordX,coordY);
      selectedStar = where((abs(starlistOrig(1,) - coordX) <= threshDetect) & \
                           (abs(starlistOrig(2,) - coordY) <= threshDetect));
      if (numberof(selectedStar) == 0) {
        gui_message,"I could not find the star ... Please re-select";
        npass -=1;
        continue;
      }
      // More than One star has been found at this location
      // Test the magnitude:
      // The brightest star should be fainter than 8 mag in R band (to be
      //     acquired by the NGS WFS)
      // The brightest star should be 2 magnitudes brighter than the others
      if (numberof(selectedStar) > 1) {
        gui_message,"More than one star correspond to this location ...";
        test_ok = 0;
        starlistInter = starlistOrig(,selectedStar);
        indexValid = where(starlistInter(5,) > -27);
        if (numberof(indexValid) < 1) {
          gui_message,"I don't have enough information to find out a suitable star, proceed to another location";
          continue;
        } else {
          if (numberof(indexValid) == 1) {
            test = where(selectedIndex == selectedStar(indexValid));
            if (numberof(test) == 0) {
              grow,selectedIndex,selectedStar(indexValid);
            } else {
              gui_message,"This star has already been selected ...proceed to another location";
              continue;
            }
            if (npass == 1) {
              starlist(,1) = starlistInter(,indexValid);
            } else grow,starlist,starlistInter(,indexValid);
            gui_message,"I only have magnitude information for one star, I selected it but proceed with caution however";
            continue;
          } else starlistInter = starlistInter(,indexValid);
        }
        magRmax = min(starlistInter(5,));
        starlistInter = starlistInter(,sort(starlistInter(5,)));

        test = where(starlistInter(5,1:-1) - 2. >= magRmax);
        if (numberof(test) > 0) {
          if (magRmax > 8.) {
            test = where(selectedIndex == selectedStar(indexValid));
            if (numberof(test) == 0) {
              grow,selectedIndex,selectedStar(indexValid);
            } else {
              gui_message,"This star has already been selected ... proceed to another location";
              continue;
            }
            if (npass == 1) {
              starlist(,1) = starlistInter(,0);
            } else grow,starlist,starlistInter(,0);
            test_ok = 1;
          } else {
            gui_message,"The brightest star is too bright for wavefront sensing";
            break;
          }
        }
        if (test_ok) {
          gui_message,"The brightest star is suitable for wavefront sensing, proceed to another location";
          continue;
        } else {
          gui_message,"I could not find a star brighter than the others ... proceed to another location";
          npass -=1;
          continue;
        }
      }
      // Only one star has been found
      // We check that it is not too bright and select it
      test = where(selectedIndex == selectedStar);
      if (numberof(test) == 0) {
        grow,selectedIndex,selectedStar;
      } else {
        gui_message,"This star has already been selected ... proceed to another location";
        continue;
      }
      if (npass == 1) {
        starlist(,1) = starlistOrig(,selectedStar);
      } else grow,starlist,starlistOrig(,selectedStar);
      gui_message,"I found the star and selected it ... proceed to another location";
    }
  } while (not_over);

  doneWithSelection = 1;
  nStarsUpdated = dimsof(starlist)(3);
  status = disp_dss_image();
  status = disp_2mass_image();

  for (nbDisp=1;nbDisp<=2;nbDisp++) {
    window,nbDisp;
    for (cptStars=1;cptStars<=nStarsUpdated;cptStars++) \
      plp,ceil(starlist(2,cptStars)/0.48)+149,ceil(starlist(1,cptStars)/0.48)+149, \
        symbol=12,width=5,size=0.8,color="red";
        // symbol=6,width=1,color="red";
  }
}

func mascot_addStarsMouse
{
  extern starlistOrig,doneWithSelection,nStarsUpdated;
  extern mascot_pixSize;

  not_over = 1;
  threshDetect = 2.;
  selectedIndex = [];

  window,1;
  gui_message,"Select stars with mouse R-click when done, M-click to remove last entry";

  do {
    res  = mouse(1,0,"");
    //pyk_status_push,"Processing...";
    c    = res(1:2);
    but  = res(10);
    if (but == 3) not_over = 0;
    if (but == 2) {
      starlist = starlist(,:-1);
      gui_message,"Last measurement taken out of star list";
      //pyk_warning,"Last measurement taken out of star list";
      continue;
    }
    if (but == 1) {
      coordX = (c(1)-149)*mascot_pixSize;
      coordY = (c(2)-149)*mascot_pixSize;
      write,coordX,coordY;
      selectedStar = where((abs(starlistOrig(1,) - coordX) <= threshDetect) & \
                           (abs(starlistOrig(2,) - coordY) <= threshDetect));
      if (numberof(selectedStar) == 0) {
        gui_message,"I could not find the star ... Please re-select";
        continue;
      }
      // More than One star has been found at this location
      // Test the magnitude:
      // The brightest star should be fainter than 8 mag in R band (to be
      //     acquired by the NGS WFS)
      // The brightest star should be 2 magnitudes brighter than the others
      if (numberof(selectedStar) > 1) {
        gui_message,"More than one star correspond to this location ...";
        test_ok = 0;
        starlistInter = starlistOrig(,selectedStar);
        indexValid = where(starlistInter(5,) > -27);
        if (numberof(indexValid) < 1) {
          gui_message,"I don't have enough information to find out a suitable star, proceed to another location";
          continue;
        } else {
          if (numberof(indexValid) == 1) {
            test = where(selectedIndex == selectedStar(indexValid));
            if (numberof(test) == 0) {
              grow,selectedIndex,selectedStar(indexValid);
            } else {
              gui_message,"This star has already been selected ... proceed to another location";
              continue;
            }
            grow,starlist,starlistInter(,indexValid);
            gui_message,"I only have magnitude information for one star, I selected it but proceed with caution";
            continue;
          } else starlistInter = starlistInter(,indexValid);
        }
        magRmax = min(starlistInter(5,));
        starlistInter = starlistInter(,sort(starlistInter(5,)));

        test = where(starlistInter(5,1:-1) - 2. >= magRmax);
        if (numberof(test) > 0) {
          if (magRmax > 8.) {
            test = where(selectedIndex == selectedStar(indexValid));
            if (numberof(test) == 0) {
              grow,selectedIndex,selectedStar(indexValid);
            } else {
              gui_message,"This star has already been selected ... proceed to another location";
              continue;
            }
            grow,starlist,starlistInter(,0);
            test_ok = 1;
          } else {
            gui_message,"The brightest star is too bright for wavefront sensing";
            break;
          }
        }
        if (test_ok) {
          gui_message,"The brightest star is suitable for wavefront sensing, proceed to another location";
          continue;
        } else {
          gui_message,"I could not find a star brighter than the others ... proceed to another location";
          continue;
        }
      }
      // Only one star has been found
      // We check that it is not too bright and select it
      test = where(selectedIndex == selectedStar);
      if (numberof(test) == 0) {
        grow,selectedIndex,selectedStar;
      } else {
        gui_message,"This star has already been selected ... proceed to another location";
        continue;
      }
      grow,starlist,starlistOrig(,selectedStar);
      gui_message,"I found the star and selected it ... proceed to another location";
    }
  } while (not_over);

  doneWithSelection = 1;
  nStarsUpdated = dimsof(starlist)(3);
  status = disp_dss_image();
  status = disp_2mass_image();

  for (nbDisp=1;nbDisp<=2;nbDisp++) {
    window,nbDisp;
    for (cptStars=1;cptStars<=nStarsUpdated;cptStars++) \
      plp,ceil(starlist(2,cptStars)/0.48)+149,ceil(starlist(1,cptStars)/0.48)+149, \
        symbol=12,width=5,size=0.8,color="red";
        // symbol=6,width=1,color="red";
  }
}


func mascot_find_best_asterism(void)
{
  extern starlist;
  extern sall;

  if (starlist==[]) {
    pyk,"glade.get_widget('bestAsterismButton').set_sensitive(1)";
    after,0.0,process_multi;
    exit,"No stars: nothing to do";
  }
  
  ns = dimsof(starlist)(0);

  if (ns>nstar_limit) {
    pyk,"glade.get_widget('bestAsterismButton').set_sensitive(1)";
    write,format="# stars (%d) too large (try <= 10)\n",ns;
    exit;
  }

  starbuf = starlist;

  sall = [];

  nast = 0;
  
  // loop on all asterism:
  
  if (ns>=3) {
    for (n1=1;n1<=ns-2;n1++) {
      for (n2=n1+1;n2<=ns-1;n2++) {
        for (n3=n2+1;n3<=ns;n3++) {
          nast++;
          starlist = starbuf(,[n1,n2,n3]);
          
          write,format="\nAsterism #%d, [%.1f,%.1f], [%.1f,%.1f], [%.1f,%.1f]\n",
            nast,starlist(1,1),starlist(2,1),starlist(1,2),
            starlist(2,2),starlist(1,3),starlist(2,3);

          // various checks on asterism:
          status = does_it_fit(starlist);
          if (status==0) {
            write,format="%s\n","Skipped. Does not fit.";
            continue;
          }
        
          sdata = mascot_compute_strehl();
          grow,sall,sdata;
          window,3;
          disp_strehl_map,sdata;

        }
      }
    }
  }

  // then still check if an asterism with 1 or 2 star lead to better perf:
  // check 2 stars asterism:
  if (ns>=2) {
    for (n1=1;n1<=ns-1;n1++) {
      for (n2=n1+1;n2<=ns;n2++) {
        nast++;
        starlist = starbuf(,[n1,n2]);
        
        write,format="\nAsterism #%d, [%.1f,%.1f], [%.1f,%.1f]\n",
          nast,starlist(1,1),starlist(2,1),starlist(1,2),starlist(2,2);
        
        // various checks on asterism:
        status = does_it_fit(starlist);
        if (status==0) {
          write,format="%s\n","Skipped. Does not fit.";
          continue;
        }
        
        sdata = mascot_compute_strehl();
        grow,sall,sdata;
        window,3;
        disp_strehl_map,sdata;
      }
    }
  }
    
  // check 1 star asterism:
  for (n1=1;n1<=ns;n1++) {
    nast++;
    starlist = starbuf(,[n1]);
      
    write,format="\nAsterism #%d, [%.1f,%.1f]\n",nast,starlist(1),starlist(2);
        
    sdata = mascot_compute_strehl();
    grow,sall,sdata;
    window,3;
    disp_strehl_map,sdata;
  }
    

  // restore value of star table
  starlist = starbuf;
  status = sort_best_asterisms();
  pyk,"glade.get_widget('bestAsterismButton').set_sensitive(1)";
  pyk,"glade.get_widget('save_current_asterism').set_sensitive(1)";

  after,0.0,process_multi;
}


func sort_best_asterisms(void)
{
  extern sall;
  extern current_ast;

  if (sall==[]) return;
  
  crit = sall.avgstrehl * (1-avg_rms_criteria) -
    sall.rmsstrehl * avg_rms_criteria;

  if (numberof(sall)>1) {
    w = sort(crit)(::-1);
    sall = sall(w);
  }
  
  current_ast = 0;
  disp_prev_next_asterism,1;
  
}

/* Main ... */
func mascot(void)
{

  if (!_pyk_proc) {
    write,"Launching mascot GUI";
    _pyk_proc = spawn(pyk_cmd, _pyk_callback);
    // gui_message,"Mascot ready !";
  } else {
    // there's already a GUI around. hence we're not going to receive
    // a signal from python to bring up windows and display. we have
    // to init the display here:
    write,"Mascot GUI already active";
  }

}

prevxy = [0,0];
mascot_ncolors=240;
pldefault,marks=0,legends=0,maxcolors=mascot_ncolors;
stop_zoom=0;
active_window = 1;
mascot_itt = 1;
image_survey_vis = dist(300);
image_2mass = dist(300);
mascot_dims = dimsof(dist(300));
mascot_cat = 0;
coord_name = "";
zoom_window = 1;
doneWithSelection = 0;
//FIXME
mascot_pixSize = 0.48;

tcsh1 = spawn("/bin/sh",callback1);
tcsh2 = spawn("/bin/sh",callback2);
tcsh3 = spawn("/bin/sh",callback3);
tcsh4 = spawn("/bin/sh",callback4);

arg = get_argv();
mascot_context = "called_from_session";
if (anyof(arg=="mascot.i")) mascot_context="called_from_shell";

// determining paths:
cds_bin_path = get_env("CDS_BIN_PATH");
if (!cds_bin_path) {
  write,format="\n%s\n"," WARNING: You haven't set up a path for the CDS binaries (aclient,";
  write,format="%s\n\n","          findnomad1, vizquery); I will assume they are in your path.";
 } else {
  if (strpart(cds_bin_path,0:0)!="/") cds_bin_path += "/";
 }

//--------------------------------
// look for python and glade files
Y_PYTHON = get_env("Y_PYTHON");
Y_GLADE  = get_env("Y_GLADE");
Y_CONF   = get_env("Y_CONF");

y_user = streplace(Y_USER,strfind("~",Y_USER),get_env("HOME"))

if (noneof(Y_PYTHON)) \
  Y_PYTHON="./:"+y_user+":"+pathform(_(y_user,Y_SITES,Y_SITE)+"python/");
if (noneof(Y_GLADE)) \
  Y_GLADE="./:"+y_user+":"+pathform(_(y_user,Y_SITES,Y_SITE)+"glade/");

// try to find mascot.py
path2py = find_in_path("mascot.py",takefirst=1,path=Y_PYTHON);
if (is_void(path2py)) {
  // not found. bust out
  pyk_error,swrite(format="Can't find mascot.py in %s.\n",Y_PYTHON);
  if (mascot_context=="called_from_shell") quit;
  error,swrite(format="Can't find mascot.py in %s.\n",Y_PYTHON);
 }
path2py = dirname(path2py);
write,format=" Found mascot.py in %s\n",path2py;

// try to find spydr.glade
path2glade = find_in_path("mascot.glade",takefirst=1,path=Y_GLADE);
if (is_void(path2glade)) {
  // not found. bust out
  pyk_error,swrite(format="Can't find mascot.glade in %s\n",Y_GLADE);
  if (mascot_context=="called_from_shell") quit;
  error,swrite(format="Can't find mascot.glade in %s\n",Y_GLADE);
 }
path2glade = dirname(path2glade);
write,format=" Found mascot.glade in %s\n",path2glade;

// spawned gtk interface
python_exec = path2py+"/mascot.py";
pyk_cmd=[python_exec,path2glade];

if (numberof(arg)>=3) mascot;
