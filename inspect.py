#!/usr/bin/env python
"""
GUI tool for inspecting grism spectrum extractions

$Date$
$Rev$

Nov 11: For team inspections:

    
Instructions for working:
- Put gui_3dhst.py in the directory where you plan to work. This is where the output
files will be created. I recommend TEST_SPECTRA_v4.1/

- Open gui_3dhst.py in your favorite text editor.

- Go to line XXX below and edit the paths for RELEASE and RGB_PATH to the correct
location on your disk. Save.

- Open a teminal window. 

- If running Ureka, type "ur_setup"

- Navigate to gui_3dhst.py directory.

Open python (type python or ipython in your terminal).
In Python:

    import gui_3dhst

	### The default behavior will pop out the ds9 window:
	gui_3dhst.go_pointing(pointing = 'goodss-05')
	
	### To avoid opening the ds9 window:
	gui_3dhst.go_pointing(pointing = 'goodss-05',show_ds9=False)
	
	     
"""
import sys
import shutil
import os
import glob

import Tkinter as tk
from PIL import ImageTk, Image

import numpy as np
import matplotlib.path as mplPath
#import pyfits
import astropy.io.fits as pyfits

noNewLine = '\x1b[1A\x1b[1M'

def examples():
    
    import glob
    #import unicorn.inspect
    import pysao
    
    os.chdir(os.getenv('THREEDHST') + '/3DHST_VariableBackgrounds/GS25')
    
    RGB_PATH = os.getenv('RELEASE') + 'v4.0/RGB/All/'
    
    ds9 = pysao.ds9()
    
    #### Objects with fits
    x = unicorn.inspect.ImageClassifier(glob.glob('ZFIT/PNG/*zfit.png'), FITS_PATH='2D/FITS/', logfile='inspect_zfit', RGB_PATH=RGB_PATH, ds9=ds9)
    
    #### Raw 2D files, including objects not fit
    x = unicorn.inspect.ImageClassifier(glob.glob('2D/PNG/*png'), FITS_PATH='2D/FITS/', logfile='inspect_raw2d', RGB_PATH=RGB_PATH, ds9=ds9)
    
    
class TextInput():
    """
    Popup for entering a comment
    """
    def __init__(self, input='---'):
        master = tk.Toplevel()
        master.geometry('250x30') # This work fine
        self.master = master
    
        self.textbox = tk.Entry(self.master)
        self.textbox.delete(0, tk.END)
        self.textbox.insert(0, input)
        
        #self.button_quit = tk.Button(self.master, text = 'Quit', command = self.finish)
        
        self.textbox.grid_rowconfigure(0, minsize=240)
        self.textbox.grid(row=0, column=0, columnspan=5)
        #self.button_quit.grid(row=0, column=5)

        self.master.bind("<Return>", self.finish)

        self.master.mainloop()

    def finish(self, event):
        #self.frame.quit()
        self.text_value = self.textbox.get()
        self.master.destroy()
        self.master.quit()
    
def go_pointing(pointing='goodss-25', RELEASE='/Volumes/Voyager/TEST_SPECTRA_v4.1/GOODS-S/INTERLACE_v4.1/', RGB_PATH='/Volumes/Voyager/TEST_SPECTRA_v4.1/RGB_goodss/',show_ds9=True):
    """
    
    Classify a 3D-HST pointing
        
    """
    import glob
    #from unicorn import inspect

    field = '-'.join(pointing.split('-')[:-1])
    
    ### Example paths for v4.0 and earlier releases
    #PNG_PATH = '%s/%s/%s-WFC3_v4.0_SPECTRA/%s/ZFIT/PNG/' %(RELEASE, field, field, pointing)
    #FITS_PATH = '%s/%s/%s-WFC3_v4.0_SPECTRA/%s/2D/FITS/' %(RELEASE, field, field, pointing)

    #PNG_PATH = '%s/%s/WFC3/%s/ZFIT/PNG/' %(RELEASE, field, pointing)
    #FITS_PATH = '%s/%s/WFC3/%s/2D/FITS/' %(RELEASE, field, pointing)
    
    PNG_PATH = FITS_PATH = RELEASE
    if show_ds9:
    	import pysao
    	ds9=pysao.ds9()
    	x = ImageClassifier(images=glob.glob(PNG_PATH+pointing+'*new_zfit.png'), RGB_PATH=RGB_PATH, FITS_PATH=FITS_PATH, logfile='%s_inspect.info' %(pointing), ds9=ds9)
    else:
    	x = ImageClassifier(images=glob.glob(PNG_PATH+pointing+'*new_zfit.png'), RGB_PATH=RGB_PATH, FITS_PATH=FITS_PATH, logfile='%s_inspect.info' %(pointing))
    return x
    
class myCheckbutton(object):
    """
    Container for CheckButton self + logging variables
    """
    def __init__(self, gui, text="(t)ilt", param=[0], hotkey='xxx'):
        self.gui = gui
        self.param = param
        self.hotkey = hotkey

        self.var = tk.IntVar()
         
        self.thebutton = tk.Checkbutton(self.gui.frame, text=text, variable=self.var, command=self.set_var)
        self.var.set(self.param[0])
        
    def set_var(self):
        self.param[self.gui.i] = self.var.get()
        #print self.gui.i
      
class mySlider(object):
    """
    Container for Scale slider self + logging variables
    """
    def __init__(self, gui, text="(a)bs./break", param=[0], hotkey='xxx', to=3):
        self.gui = gui
        self.param = param
        self.hotkey = hotkey
        self.to = to        
        
        self.var = tk.IntVar()
        
        self.theslider = tk.Scale(self.gui.frame, to=to, resolution=1, orient=tk.HORIZONTAL, variable=self.var, label=text)
        self.var.set(param[0])
        
    def set_var(self):
        self.param[self.gui.i] = self.var.get()
    
    
class ImageClassifier():
    """
    Main classifier tool for 3D-HST fits
    """
    def __init__(self, images=[], logfile='inspect_raw.info', load_log=True, 
                 master=None):
        """
        GUI tool for inspecting grism redshift fits
        
         x = unicorn.inspect.ImageClassifier(images=glob.glob('Specz/GOODS-S-25*zfit.png'), RGB_PATH=RGB_PATH, FITS_PATH='Specz/')
                  
         """
        if len(images) == 0:
            print 'No images specified'
            return False
            
        if not os.path.exists(images[0]):
            print 'First image not found (%s), is path correct?' %(images[0])
            return False
        
        ##### Add .fits to filename and make backup if necessary
        self.logfile = logfile
        if not self.logfile.lower().endswith('.fits'):
            self.logfile += '.fits'
        
        if os.path.exists(self.logfile):
            bk = glob.glob(self.logfile+'.backup*')
            if len(bk) > 0:
                bkup_file = self.logfile + '.backup.%03d' %(len(bk))
            else:
                bkup_file = self.logfile + '.backup'
                
            shutil.copy(self.logfile, bkup_file)
            print 'Made copy of %s -> %s' %(self.logfile, bkup_file)
        
        ####### Initialize parameters
        self.params = {}        
        self.images = images
        
        self.marked_reads = None
        self.NREAD = 14
        
        ### Polygons for reads
        x0 = y0 = 12
        px = py = 6
        dx = dy = 241
        xi = np.array([0,1,1,0])
        yi = np.array([0,0,1,1])
        
        c = 0
        self.read_polygons = []
        for j in range(4):
            for i in range(4):
                c += 1
                if c > self.NREAD:
                    break
                else:
                    polyx = x0+i*(px+dx)+xi*dx
                    polyy = y0+j*(py+dy)+yi*dy
                    poly = np.array([polyx, polyy]).T
                    self.read_polygons.append(mplPath.Path(poly))
                    
        if os.path.exists(self.logfile) & load_log:
            self.read_fits()
            
        self.N = len(self.images)

        for key in ['satellite', 'earth', 'other', 'kill', 'seen']:
            if key not in self.params.keys():
                self.params[key] = np.zeros(self.N, dtype=np.int)
        
        if self.marked_reads is None:
            self.marked_reads = np.zeros((self.N, self.NREAD), dtype=int)
        
        if 'comment' not in self.params.keys():
            self.params['comment'] = ['---' for i in range(self.N)]
                                        
        self.i = 0
        self.master = master
        self.setup_gui()
        
    def setup_gui(self):
        
        ####### Initialize GUI
        if self.master is None:
            master = tk.Toplevel()
            master.geometry('1050x1200') 
            self.master = master
         
        self.frame = tk.Frame(self.master)
        self.frame.pack()
                
        #### Image Panels
        imageFile = Image.open(self.images[0]).resize((1000,1000))
        im = ImageTk.PhotoImage(imageFile)        
        self.panel = tk.Label(self.frame , image=im, cursor='target')
        
        #### Keypress binding
        self.master.bind("<Key>", self.keypress_event)
        
        #### Mouse binding - right click
        self.master.bind("<Button-2>", self.right_click_event)
        
        ### For logging slider movements
                
        ######
        ### Navigation buttons
        ###
        self.button_quit = tk.Button(self.frame, text = '(q)uit', command = self.finish)
        self.button_log = tk.Button(self.frame, text = 'Log to (F)ile', command = self.write_fits)
        
        self.button_prev = tk.Button(self.frame, text = '(p)rev', command = self.img_prev)
        self.button_next = tk.Button(self.frame, text = '(n)ext', command = self.img_next)
        
        self.buttons = {}
            
        #### Satellite trail
        self.buttons['satellite'] = myCheckbutton(self, text='(s)atellite', param=self.params['satellite'], hotkey='s')
        
        #### Earth glow
        self.buttons['earth'] = myCheckbutton(self, text='(e)arthglow', param=self.params['earth'], hotkey='e')
        
        #### Other
        self.buttons['other'] = myCheckbutton(self, text='(o)ther', param=self.params['other'], hotkey='o')
        
        #### Kill
        self.buttons['kill'] = myCheckbutton(self, text='(k)ill', param=self.params['kill'], hotkey='k')
        
        ### Object seen already?
        self.buttons['seen'] = myCheckbutton(self, text='SEEN?', param=self.params['seen'], hotkey='xxx')

        ### Read holder
        self.rvar = tk.StringVar()
        self.r_text = tk.Label(self.frame, textvariable=self.rvar)
        #self.e_comment.configure(relief=tk.RIDGE)
        self.rvar.set("Test")        
        
        ### Comment holder
        self.tvar = tk.StringVar()
        self.e_comment = tk.Label(self.frame, textvariable=self.tvar)
        self.e_comment.configure(relief=tk.RIDGE)
        self.tvar.set(self.params['comment'][0])        
        
        #####################        
        ##### Set up the grid for the GUI elements
        self.button_next.grid(row=1, column=4, columnspan=1)
        self.button_prev.grid(row=1, column=5, columnspan=1)
        
        self.button_log.grid(row=2, column=4, columnspan=1)
        self.button_quit.grid(row=2, column=5, columnspan=1)

        self.buttons['satellite'].thebutton.grid(row=1, column=1)
        self.buttons['earth'].thebutton.grid(row=1, column=2)

        self.buttons['other'].thebutton.grid(row=2, column=1)
        self.buttons['kill'].thebutton.grid(row=2, column=2)

        self.r_text.grid(row=3, column=1, columnspan=2)
        
        self.e_comment.grid(row=3, column=4, columnspan=2)
        
        self.canvas = tk.Canvas(self.frame, background="white", height=5)
        self.line_tags = []
        self.draw_lines()
        
        self.panel.grid(row=0, column=0, columnspan=6)
        #self.canvas.grid(row=0, column=0, columnspan=6)
        self.canvas.place(relx=0, rely=0, relwidth=1)
        self.panel.grid()
                
        self.master.mainloop()
    
    def right_click_event(self, event):
        """
        Mouse was clicked on a stack spectrum figure
        """
        #print "Clicked at %.1f %.1f" %(event.x, event.y)
        #oval = self.canvas.create_oval(event.x-5, 10, event.x+5, 0, fill="red", outline="blue", width=1, tags="line tag")
        #box = self.canvas.create_rectangle(event.x-5, -10, event.x+5, 10, fill="red", outline="white", width=1, tags="line tag")
        
        #self.canvas.tag_raise(self.panel)
        for j in range(self.NREAD):
            if self.read_polygons[j].contains_point((event.x, event.y)):
                print 'Read #%d!' %(j+1)
                
                self.marked_reads[self.i, j] = (not self.marked_reads[self.i, j])*1
                self.draw_lines()
                return None
                
            # if np.abs(self.marked_reads[self.i, j] - event.x) < 8:
            #     self.marked_reads[self.i, j] = 0.
            #     self.draw_lines()
            #     return None
            
        # j = 0
        # while (self.marked_reads[self.i, j] > 0) & (j < 14):
        #     j += 1
        # 
        # self.marked_reads[self.i, j] = event.x*1.
        # 
        # #print j, self.marked_reads[self.i,:]
        # self.draw_lines()
        return None
        
    def draw_lines(self):
        
        # self.canvas.delete("line tag")
        # for j in range(len(self.line_tags)):
        #     p = self.line_tags.pop()
        #     #print 'kill: ', p
        #     self.canvas.delete(p)
        
        marked = []    
        for j in range(self.NREAD):
            x = self.marked_reads[self.i, j]
            if x > 0:
                # verts = self.read_polygons[j].vertices
                # print j, verts
                # tag = self.canvas.create_rectangle(verts[0][0], verts[0][1], verts[2][0], verts[2][1], fill="red", outline="green", width=2, tags="line tag")
                # #print tag
                # self.line_tags.append(tag)
                marked.append(j+1)
        
        marked.sort()
        self.rvar.set('Marked: %s' %(marked))
        
        #print 'Raise!'
        #self.canvas.tag_raise("line tag")
        #self.canvas.tag_lower(self.panel)
        
    def keypress_event(self, event):
        key = event.char
        #print 'Keyboard: "%s"' %(key)
        #self.line_flags[self.i] = (not self.params['line'][self.i])*1
        
        #### Buttons
        for id in self.buttons.keys():
            if key == self.buttons[id].hotkey:
                button = self.buttons[id]
                button.param[self.i] = (not button.param[self.i])*1
                button.var.set(button.param[self.i])
                #print button.param[self.i]
                return True
        
        # #### Sliders
        # for id in self.sliders.keys():
        #     if key == self.sliders[id].hotkey:
        #         slider = self.sliders[id]
        #         slider.param[self.i] = ((slider.param[self.i]+1) % (slider.to+1))
        #         slider.var.set(slider.param[self.i])
        #         #print slider.param[self.i]
        #         return True

        #### Other keys
        if key == 'c':
            #comment = raw_input('Comment: ')
            ctext = TextInput(self.params['comment'][self.i])
            self.tvar.set(ctext.text_value)
            self.params['comment'][self.i] = ctext.text_value
            #print comment
            
        elif key == 'n':
            ### next
            self.img_next()
        
        elif key == 'p':
            ### previous
            self.img_prev()
        
        elif key == 'f':
            ### Force write the output file
            self.write_fits()
        
        elif key == 'q':
            ### quit
            self.finish()

        elif key == 'N':
            ### Go to next unseen
            self.set_seen()
            self.params['comment'][self.i] = self.tvar.get()
            
            while (self.i < self.N-1):
                #print self.i
                if self.params['seen'][self.i] == 1:
                    self.i += 1
                else:
                    break

            #self.img_next()
            self.load_image()

        elif key == '0':
            ### Go to beginning
            self.set_seen()
            self.params['comment'][self.i] = self.tvar.get()
            
            self.i = 0
            self.load_image()
                    
        elif key == '!':
            ### Go until find '!' in a comment
            self.set_seen()
            self.params['comment'][self.i] = self.tvar.get()

            ### At end
            if self.i > self.N-2:
                return False
            
            self.i += 1
            
            while (self.i < self.N-1):
                print self.i
                if (self.params['seen'][self.i] == 1) & ('!' not in self.params['comment'][self.i]):
                    self.i += 1
                else:
                    break
            
            self.load_image()
            
            # while (self.params['seen'][self.i] == 1) & ('!' not in self.params['comment'][self.i]):
            #     self.img_next()
        
        elif (key in 'KOES') & (key is not ''):
            ### Go to next non-zero flag for different keys
            key_param = {'K':self.params['kill'], 'O':self.params['other'], 'E':self.params['earth'], 'S':self.params['satellite']}
            self.next_nonzero(key_param[key])
            
        elif key == '?':
            print """
    Additional keys:

      'c':  Open comment box, type <tab> to edit and <enter> when done.
      '9':  Open 2D FITS extensions in DS9, if pysao available.
      'N':  skip through through all 'seen'
      '0':  go back to first
      '!':  go to next object with '!' in the comment
      'y':  go to next unambigous
      '?':  this message
      
"""
        else:
            print 'Hotkey (%s) not bound.' %(key)
    
    def next_nonzero(self, param=None):
        """
        Go to next object where self.param is nonzero
        """
        self.set_seen()
        self.params['comment'][self.i] = self.tvar.get()

        ### At end
        if self.i > self.N-2:
            return False
        
        self.i += 1
        
        while (self.i < self.N-1):
            print self.i
            if (self.params['seen'][self.i] == 1) & (param[self.i] == 0):
                self.i += 1
            else:
                break
        
        self.load_image()
            
    def set_seen(self):
        #print 'Seen flag: %d' %(self.svar.get())
        self.params['seen'][self.i] = 1
    
    def img_next(self):
        self.set_seen()
        self.params['comment'][self.i] = self.tvar.get()
        
        if self.i == self.N-1:
            print 'Already at last image'
            return False
                    
        self.i += 1
        self.load_image()
        self.draw_lines()
                    
        return True
    #
    def img_prev(self):
        self.set_seen()
        self.params['comment'][self.i] = self.tvar.get() #self.e_comment.get()

        if self.i == 0:
            print 'Already at first image'
            return False
                    
        self.i -= 1
        self.load_image()
        self.draw_lines()
        
        return True
        
    def load_image(self):
        
        #print '%d  %s' %(self.i, self.images[self.i])
        #print '%d %d %d' %(self.params['line'][self.i], self.params['unamb'][self.i], self.params['star'][self.i])
        
        ramp_file = self.images[self.i].replace('ramp.png','ramp.dat')
        pop_reads = check_background_SN(ramp_file=ramp_file, show=False)
        
        print '%s: %d of %d %s' %(self.images[self.i], self.i+1, self.N, pop_reads)
        
        for key in self.buttons.keys():
            button = self.buttons[key]
            button.var.set(button.param[self.i])
        
        # for key in self.sliders.keys():
        #     slider = self.sliders[key]
        #     slider.var.set(slider.param[self.i])
        
        #self.e_comment.delete(0, tk.END)
        #self.e_comment.insert(0, self.params['comment'][self.i])
        self.tvar.set(self.params['comment'][self.i])
        
        imageFile = Image.open(self.images[self.i]).resize((1000,1000))
        im = ImageTk.PhotoImage(imageFile)
        self.panel.configure(image = im)
        self.panel.image = im
        
        #self.panel_rgb = tk.Label(self.frame , image=im_rgb)
        
    def write_fits(self):
        """
        Write the FITS log
        """
        
        import time
        import getpass
        
        formats = {}
        formats['bool'] = 'L'
        formats['int16'] = 'I'
        formats['int32'] = 'J'
        formats['int64'] = 'K'
        formats['float32'] = 'E'
        formats['float64'] = 'D'
        
        formats['>i8'] = 'K'
        formats['>f8'] = 'D'
        
        #### Make the table columns, translating numpy data types to "TFORM"
        coldefs = []
        TFORM = 'A'+str(np.array(self.images).dtype).split('S')[1]
        
        coldefs.append(pyfits.Column(name='images', array=np.array(self.images), format=TFORM))
        
        for column in self.params.keys():
            if column == 'comment':
                coldata = np.array(self.params['comment'])
            else:
                coldata = self.params[column]
            #
            dtype = str(coldata.dtype)
            #print column, dtype
            if dtype in formats.keys():
                TFORM=formats[dtype]
            else:
                if 'S' not in dtype:
                    print 'Unrecognized data type in: %s' %(dtype)
                    return False
                #
                TFORM = 'A'+dtype.split('S')[1]
            #
            #data = self.params[column]
            if '>' in dtype:
                cast_types = {'>i8':np.int64, '>f8':np.float64}
                coldata = np.cast[cast_types[dtype]](coldata)
            #
            coldefs.append(pyfits.Column(name=column, array=coldata, format=TFORM))
        
        #### Done, now make the binary table
        tbhdu = pyfits.BinTableHDU().from_columns(coldefs)
        
        linehdu = pyfits.ImageHDU(data=self.marked_reads, name='FLAGGED')
        
        #### Primary HDU
        hdu = pyfits.PrimaryHDU()
        thdulist = pyfits.HDUList([hdu, tbhdu, linehdu])

        #### Add modification time of "infile" to FITS header
        infile_mod_time = time.strftime("%m/%d/%Y %I:%M:%S %p",
                            time.localtime()) # os.path.getmtime(self.filename)))
        
        thdulist[0].header['MODTIME'] =  infile_mod_time
        thdulist[0].header['USER'] = getpass.getuser()
            
        thdulist.writeto(self.logfile, clobber=True)
        
        print 'Log to file %s' %(self.logfile)
        
    def read_fits(self):
        """
        Read already saved output
        """
        
        im = pyfits.open(self.logfile)
        if len(im) == 3:
            try:
                self.marked_reads = im['FLAGGED'].data
            except:
                self.marked_reads = None
        
        #
        tab = im[1].data
        print "Read log %s from %s (%s)" %(self.logfile, im[0].header['USER'], im[0].header['MODTIME'])
        
        colnames = tab.columns.names
        
        try:
            #### If images in the input list aren't in the file that 
            #### was read, add them to the table
            from astropy.table import Table
            tab = Table.read(self.logfile)
            colnames = tab.columns
            read_images = tab['images']
            if self.marked_reads is None:
                self.marked_reads = np.zeros((len(tab), self.NREAD), dtype=int)
                
            Nadd = 0
            for image in self.images:
                if image not in read_images:
                    Nadd += 1
                    tab.add_row()
                    tab['images'][-1] = image
                    tab['comment'][-1] = '---'
            
            if Nadd > 0:
                self.marked_reads = np.append(self.marked_reads, np.zeros((Nadd, self.NREAD), dtype=int), axis=0)
                
        except:
            #### No astropy?
            print 'No astropy found.  Forcing image list from %s.' %(self.logfile)
            pass
            
        self.images = tab['images']
            
        for c in colnames:
            #print c
            if c == 'images':
                continue

            if c == 'comment':
                self.params[c] = list(tab[c])
                continue
            #    
            self.params[c] = tab[c]
            
    def finish(self):
        #self.frame.quit()
        #self.write_log()
        self.write_fits()
        self.master.destroy()
        self.master.quit()
#
def check_background_SN(ramp_file='id1ketvxq_ramp.dat', show=False):
    """
    Check how S/N evolves with adding reads with increasing background
    """
    from astropy.table import Table
    import matplotlib.pyplot as plt
    import numpy as np
    
    time, bg = np.loadtxt(ramp_file, unpack=True)
    t0 = np.diff(time)[0]
    time, bg = time[1:], bg[1:] # skip first 2.9 s read

    if len(bg) <= 2:
        return []
        
    dt = np.append(t0, np.diff(time))
    
    s = 1. # test count rate
    so = np.argsort(bg)
    
    NREAD = len(bg)
    reads = np.arange(NREAD)+1
    
    fluence = np.cumsum((s*dt)[so])
    rms = np.sqrt(np.cumsum((bg*dt)[so]))
    rms_optimal = np.sqrt(np.cumsum((bg[so][2]*dt)[so]))

    max_ix = np.argmax(fluence/rms)
    if max_ix < (NREAD-1):
        pop_reads = list(reads[so][max_ix+1:])
    else:
        pop_reads = []
    
    ### where observed rms < expected RMS for flat background
    other_bad = list(reads[so][rms/rms_optimal > 1.25])
    pop_reads = np.cast[int](np.unique(np.hstack((pop_reads, other_bad))))
    pop_reads = list(pop_reads)
    
    if len(pop_reads) > 0:
        # np.savetxt(ramp_file.replace('ramp.dat', 'ramp.pop.png'), pop_reads,
        #            fmt='%d')
        fp = open('/tmp/'+ramp_file.replace('ramp.dat', 'ramp.pop.dat'),'w')
        for r in pop_reads:
            fp.write('%d\n' %(r))
        fp.close()
        
    if show:
        plt.ioff()
        
        fig = plt.figure(figsize=[6,3])
        ax = fig.add_subplot(111)
        ax.plot(time, bg, color='0.8', linewidth=6, zorder=1)
        si = 60
        ax.scatter(time, bg, color='w', s=si, zorder=8)
        ax.scatter(time, bg, color='k', s=0.5*si, zorder=9)
        if len(pop_reads) > 0:
            ix = np.array(pop_reads)
            ax.scatter(time[ix-1], bg[ix-1], color='r', zorder=10, s=0.5*si)
        
        ax.grid()
        ax.set_title(ramp_file)
        ax.set_xlabel('time'), ax.set_ylabel('background')
        
        fig.tight_layout(pad=0.1)
        
        fig.savefig('/tmp/'+ramp_file.replace('.dat', '.pop.png'))
        plt.close()
        
    return pop_reads
        
if __name__ == "__main__":
    files = sys.argv[1:]
    widget = ImageClassifier(images=files)
