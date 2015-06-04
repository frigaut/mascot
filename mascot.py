#!/usr/bin/env python2
# mascot.py

import gtk
import gtk.glade
import sys
import gobject
import os, fcntl, errno

class mascot:
   
   def destroy(self, wdg, data=None):
      self.py2yo('quit')
      raise SystemExit
#      gtk.main_quit()
      
   def __init__(self,mascottop):
      self.mascottop = mascottop
      self.usercmd = 'STOP'

      self.glade = gtk.glade.XML(os.path.join(self.mascottop,'mascot.glade')) 
      self.window = self.glade.get_widget('window1')
      self.statusbar = self.glade.get_widget('statusbar')
      self.progressbar = self.glade.get_widget('progressbar')
          
     # handle destroy event
      if (self.window):
         self.window.connect('destroy', self.destroy)
         
      # self.glade.signal_autoconnect(dic)
      self.glade.signal_autoconnect(self)
      
      # set stdin non blocking, this will prevent readline to block
      fd = sys.stdin.fileno()
      flags = fcntl.fcntl(fd, fcntl.F_GETFL)
      fcntl.fcntl(fd, fcntl.F_SETFL, flags | os.O_NONBLOCK)
      
      # add stdin to the event loop (yorick input pipe by spawn)
      gobject.io_add_watch(sys.stdin,gobject.IO_IN|gobject.IO_HUP,self.yo2py,None)

      # update parameters from yorick:
      #self.py2yo('gui_update')

      #self.glade.get_widget('wfs_and_dms').hide()
      #ebox = self.glade.get_widget('eventbox1')
      #ebox.connect('key-press-event',self.on_eventbox1_key_press)
      
      self.glade.get_widget('CatalogSelect1').set_active(0)
      self.glade.get_widget('comboboxentry1').set_active(0)
      
      # run
      gtk.main()

   doing_zoom=0
   done_init=0
   info_type=1
   n_stars=3
      
   def on_skin1_activate(self, *args):
      rc_defs = list(gtk.rc_get_default_files())
      # rc_xtra = os.path.join(mysttop,'resources/mcao-eng.gtkrc')
      rc_xtra = '/export/home/mcao/mcao-mascot/resources/mcao-eng.gtkrc'
      if rc_xtra in rc_defs: rc_defs.remove(rc_xtra)
      gtk.rc_set_default_files(rc_defs)
      settings = gtk.settings_get_default()
      gtk.rc_reparse_all_for_settings(settings, force_load=True)

   def on_skin2_activate(self, *args):
      rc_defs = list(gtk.rc_get_default_files())
      # rc_xtra = os.path.join(mysttop,'resources/mcao-eng.gtkrc')
      rc_xtra = '/export/home/mcao/mcao-mascot/resources/mcao-eng.gtkrc'
      if rc_xtra not in rc_defs: rc_defs.append(rc_xtra)
      gtk.rc_set_default_files(rc_defs)
      settings = gtk.settings_get_default()
      gtk.rc_reparse_all_for_settings(settings, force_load=True)
      

   def on_shortcut_help_activate(self,wdg):
      self.py2yo('mascot_shortcut_help')

   def on_about_activate(self,wdg):
     dialog = self.glade.get_widget('aboutdialog')
     dialog.run()
     dialog.hide()

   def on_clear_cache_activate(self,wdg):
      self.py2yo('clear_cache')


   #
   # Yorick to Python Wrapper Functions
   #
      
   def on_quit_activate(self,*args):
      self.py2yo('quit')
      raise SystemExit

   def on_window1_map_event(self,wdg,*args):
      drawingarea1 = self.glade.get_widget('dssdrawing')
      mwid1 = drawingarea1.window.xid;
      drawingarea2 = self.glade.get_widget('2massdrawing')
      mwid2 = drawingarea2.window.xid;
      drawingarea3 = self.glade.get_widget('strehldrawing')
      mwid3 = drawingarea3.window.xid;
      drawingarea4 = self.glade.get_widget('zoomdrawing')
      # mwid4 = drawingarea4.window.xid;
      # self.py2yo('mascot_win_init %d %d %d %d' % (mwid1,mwid2,mwid3,mwid4))
      self.py2yo('mascot_win_init %d %d %d %d' % (mwid1,mwid2,mwid3,0))

   def on_ebox_dss_enter_notify_event(self,wdg,*arg):
      self.window.set_focus(wdg)
      
   def on_ebox_dss_key_press_event(self,wdg,event):
      # sys.stderr.write('d pressed\n')
      if (event.string=='d'):
         self.py2yo('delete_star_under_cursor')
      if (event.string=='i'):
         self.py2yo('info_star_under_cursor')
      if (event.string=='c'):
         self.py2yo('recenter_here')
      if (event.string=='o'):
         self.py2yo('recenter_here')
      if (event.string=='?'):
         self.py2yo('mascot_shortcut_help')
      if (event.string=='/'):
         self.py2yo('mascot_shortcut_help')
      
   # def on_dssdrawing_enter_notify_event(self,wdg,*args):
   #    self.glade.get_widget('eventbox1').grab_focus()
   #    if (self.doing_zoom==0):
   #       self.py2yo('pyk_set stop_zoom 0')
   #       self.py2yo('pyk_set zoom_window 1')
   #       self.py2yo('mascot_zoom')
   #       self.doing_zoom=1

   # def on_dssdrawing_leave_notify_event(self,wdg,*args):
   #    self.doing_zoom=0
   #    self.py2yo('pyk_set stop_zoom 1')
            
   # def on_2massdrawing_enter_notify_event(self,wdg,*args):
   #    self.glade.get_widget('eventbox2').grab_focus()
   #    if (self.doing_zoom==0):
   #       self.py2yo('pyk_set stop_zoom 0')
   #       self.py2yo('pyk_set zoom_window 2')
   #       self.py2yo('mascot_zoom')
   #       self.doing_zoom=1
         
   # def on_2massdrawing_leave_notify_event(self,wdg,*args):
   #    self.doing_zoom=0
   #    self.py2yo('pyk_set stop_zoom 1')
      
   # def on_strehldrawing_enter_notify_event(self,wdg,*args):
   #    self.glade.get_widget('strehldrawing').grab_focus()
   #    if (self.doing_zoom==0):
   #       self.py2yo('pyk_set zoom_window 3')
   #       self.py2yo('mascot_zoom')
   #       self.doing_zoom=1

   # def on_strehldrawing_leave_notify_event(self,wdg,*args):
   #    self.doing_zoom=0
   #    self.py2yo('pyk_set stop_zoom 1')
      
   def on_coordButton1_pressed(self,wdg,*args):
      self.info_type=2

   def on_targetButton1_pressed(self,wdg,*args):
      self.info_type=1


   def on_targetCat_changed(self,wdg):
      targetCat = wdg.get_active_text()
      if (targetCat=="NGC/IC Catalog"):
         self.py2yo('pyk_set mascot_cat 1')
      elif (targetCat=="Globular Clusters Catalog"):
         self.py2yo('pyk_set mascot_cat 2')
      elif (targetCat=="Galaxies Catalog"):
         self.py2yo('pyk_set mascot_cat 3')

   def on_starsCat_changed(self,wdg):
      targetCat = wdg.get_active_text()
      if (targetCat=="Nomad1"):
         self.py2yo('pyk_set mascot_starcat 1')
      elif (targetCat=="Guide Stars Catalog"):
         self.py2yo('pyk_set mascot_starcat 2')

   def on_SearchTargetButton_pressed(self,wdg):
      target = self.glade.get_widget('targetEntry').get_text()
      if len(target):
         self.window.window.set_title("mascot: %s" % target);
         # self.py2yo('pyk_set coord_name "%s"' % target)
         self.py2yo('mascot_search_target "%s"' % target)
         self.glade.get_widget('SearchTargetButton').set_sensitive(0)
         self.glade.get_widget('targetEntry').set_sensitive(0)
         self.glade.get_widget('SearchStarsButton').set_sensitive(0)
         self.glade.get_widget('star_filters').set_sensitive(0);
         self.glade.get_widget('bestAsterismButton').set_sensitive(0)
      
   def on_targetEntry_activate(self,wdg):
      target = self.glade.get_widget('targetEntry').get_text()
      if len(target):
         self.window.window.set_title("mascot: %s" % target);
         # self.py2yo('pyk_set coord_name "%s"' % target)
         self.py2yo('mascot_search_target "%s"' % target)
         self.glade.get_widget('SearchTargetButton').set_sensitive(0)
         self.glade.get_widget('targetEntry').set_sensitive(0)
         self.glade.get_widget('SearchStarsButton').set_sensitive(0)
         self.glade.get_widget('star_filters').set_sensitive(0);
         self.glade.get_widget('bestAsterismButton').set_sensitive(0)

   def on_minmag_activate(self,wdg):
      self.py2yo('pyk_set mag_min_threshold %s' % wdg.get_text())
      self.py2yo('mascot_search_stars')
      self.glade.get_widget('SearchStarsButton').set_sensitive(0)
      self.glade.get_widget('star_filters').set_sensitive(0);
      self.glade.get_widget('bestAsterismButton').set_sensitive(0)
      
   def on_maxmag_activate(self,wdg):
      self.py2yo('pyk_set mag_max_threshold %s' % wdg.get_text())
      self.py2yo('mascot_search_stars')
      self.glade.get_widget('SearchStarsButton').set_sensitive(0)
      self.glade.get_widget('star_filters').set_sensitive(0);
      self.glade.get_widget('bestAsterismButton').set_sensitive(0)
            
   def on_SearchStarsButton_pressed(self,wdg):
      self.py2yo('mascot_search_stars')
      self.glade.get_widget('SearchStarsButton').set_sensitive(0)
      self.glade.get_widget('star_filters').set_sensitive(0);
      self.glade.get_widget('bestAsterismButton').set_sensitive(0)
      
   def on_SelectStarsButton_pressed(self,wdg):
      self.py2yo('mascot_select_stars')
  
   def on_ExcludeStarsButton_pressed(self,wdg):
      self.py2yo('mascot_excludeStars %d' % (n_stars))

   def on_addStar_mouse_clicked(self,wdg):
      self.py2yo('mascot_addStarsMouse')
      
   def on_bestAsterismButton_pressed(self,wdg):
      self.py2yo('mascot_find_best_asterism')
      self.glade.get_widget('bestAsterismButton').set_sensitive(0)

#   def on_strehlButton_pressed(self,wdg):
#      self.py2yo('mascot_compute_strehl')
      
   def on_prev_ast_pressed(self,wdg):
      self.py2yo('disp_prev_next_asterism -1')
   
   def on_next_ast_pressed(self,wdg):
      self.py2yo('disp_prev_next_asterism 1')

   def on_save_current_asterism_pressed(self,wdg):
      self.py2yo('print_asterism current_ast 1')

   def on_avg_rms_crit_value_changed(self,wdg):
     self.py2yo('pyk_set avg_rms_criteria %f' % wdg.get_value())
     self.py2yo('sort_best_asterisms')

   def on_minmag_changed(self,wdg):
     self.py2yo('pyk_set mag_min_threshold %s' % wdg.get_text())
     
   def on_maxmag_changed(self,wdg):
     self.py2yo('pyk_set mag_max_threshold %s' % wdg.get_text())
     
   def on_nstar_limit_changed(self,wdg):
     self.py2yo('pyk_set nstar_limit %s' % wdg.get_text())
     
   def y_set_checkbutton(self,name,val):
      self.glade.get_widget(name).set_active(val)

   def pyk_status_push(self,id,txt):
      self.glade.get_widget('statusbar').push(id,txt)
      
   def pyk_status_pop(self,id):
      self.glade.get_widget('statusbar').pop(id)
      
   def pyk_error(self,msg):
      dialog = gtk.MessageDialog(type=gtk.MESSAGE_ERROR,buttons=gtk.BUTTONS_OK,message_format=msg)
      dialog.run()
      dialog.destroy()

   def pyk_info(self,msg):
      dialog = gtk.MessageDialog(type=gtk.MESSAGE_INFO,buttons=gtk.BUTTONS_OK,message_format=msg)
      dialog.run()
      dialog.destroy()

   def pyk_warning(self,msg):
      dialog = gtk.MessageDialog(type=gtk.MESSAGE_WARNING,buttons=gtk.BUTTONS_OK,message_format=msg)
      dialog.run()
      dialog.destroy()

   def on_invert_toggled(self,wdg):
      if (self.done_init):
         self.py2yo('pyk_set mascot_invertlut %d' % wdg.get_active())
         self.py2yo('mascot_lut')

   def on_unzoom_activate(self,wdg):
      self.py2yo('unzoom_all')
      
   def on_limits_clicked(self,wdg):
      self.py2yo('unzoom %d' % self.active_window)      
      self.py2yo('do_limits')
      
   def on_mascot_help_activate(self,wdg):
      self.py2yo('mascot_shortcut_help')

   #
   # minimal wrapper for yorick/python communication
   #
   
   def py2yo(self,msg):
      # sends string command to yorick's eval
      sys.stdout.write(msg+'\n')
      sys.stdout.flush()
   
   def yo2py(self,cb_condition,*args):
      if cb_condition == gobject.IO_HUP:
         raise SystemExit, "lost pipe to yorick"
      # handles string command from yorick
      # note: inidividual message needs to end with /n for proper ungarbling
      while 1:
         try:
            msg = sys.stdin.readline()
            msg = "self."+msg
            exec(msg)
         except IOError, e:
            if e.errno == errno.EAGAIN:
               # the pipe's empty, good
               break
            # else bomb out
            raise SystemExit, "yo2py unexpected IOError:" + str(e)
         except Exception, ee:
            raise SystemExit, "yo2py unexpected Exception:" + str(ee)
      return True

   def set_cursor_busy(self,state):
      if state:
         self.window.window.set_cursor(gtk.gdk.Cursor(gtk.gdk.WATCH))
      else:
         self.window.window.set_cursor(gtk.gdk.Cursor(gtk.gdk.LEFT_PTR))
         
if len(sys.argv) != 2:
   print 'Usage: mascot.py path_to_mascot'
   raise SystemExit

mascottop = str(sys.argv[1])
top = mascot(mascottop)
