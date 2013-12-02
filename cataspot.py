#!/bin/env python
"""

catapot: a program for collecting misc. data from a series of images

for more info, see the wiki page


djo, 7/13

"""


# ------------------------- imports -------------------------
# std lib
import collections
import getpass
import json
import optparse
import os
import time

# image and math
import Image
import ImageTk
import ImageDraw
import numpy

# GUI
import Tkinter
import Pmw
import tkFileDialog as tkFD
import tkMessageBox as tkMB



# ------------------------- constants -------------------------

__version__ = "0.5"

defaulttitle = "Cataspot"


# ----- file stuff

# ongoing save file (hardcoded name):
savefilename = "cataspot.json"
savefileversion = 1
savefiledescription = "cataspot data file"


# ----- graphics stuff
blankimage = Image.new('RGB', (1024, 1024))

# contrast: for now, set max to some fraction of data max:
maxcontrast = 0.1

# ----- tool stuff
Marker = collections.namedtuple("Marker", ["kind", "size", "color"])

toolmarkerdict = {
    "centering reflections": Marker("+", 10, (255, 0, 0)),
    "beam center": Marker("(+)", 10, (255, 0, 0)),
    "reference points": Marker("+", 10, (0, 255, 255)),
    "beam stop center": Marker("+", 10, (255, 255, 0)),
    "spots": Marker("x", 10, (0, 255, 0)),
    }
# repeat the list in the order we want:
toollist = [
    "centering reflections", 
    "reference points", 
    "beam stop center", 
    "spots", 
    ]
defaulttool = toollist[0]


# ----- misc.
# template for status bar label; fill with x, y, image value
coordinatetemplate = "%5d  %5d  %5d"


# ------------------------- class CataspotApp -------------------------
class CataspotApp(object):
    """
    the main application
    """
    # ......................... __init__ .........................
    def __init__(self, options):
        """
        
        input:
        """

        # bookkeeping
        self.options = options


        # no data initially; all results are held in this dict, which
        #   will eventually end up in json files:
        self._data = {}

        self._imagedirectory = ""


        # ------------------------- start up the gui -------------------------
        self.root = Tkinter.Tk()
        Pmw.initialise(self.root)
        self.root.title(defaulttitle)


        self.setupmenus()

        # make the window big, but fit it within the current screen
        screenx = min(self.root.winfo_screenwidth(), 1600)
        screeny = min(self.root.winfo_screenheight(), 1200)
        self.root.geometry("%sx%s+50+50" % (screenx - 200, screeny - 200))

        # big frame and status bar beneath:
        mainframe = Tkinter.Frame(self.root)
        mainframe.pack(side='top', expand='yes', fill='both')

        # status bar:
        self.coordinatelabel = Tkinter.Label(self.root, relief='sunken')
        self.coordinatelabel.configure(text=coordinatetemplate % (0, 0, 0))
        self.coordinatelabel.pack(side='top', fill='x', anchor='w', padx=4, pady=4)





        # ----- file list area
        fileframe = Tkinter.Frame(mainframe)
        fileframe.pack(side='left', fill='both')

        self.folderentry = Pmw.EntryField(fileframe,
            labelpos='n',
            label_text="Directory",
            )
        self.folderentry.pack(side='top', padx=5, pady=5, fill='x')

        self.pendinglistbox = Pmw.ScrolledListBox(fileframe,
            labelpos='n',
            label_text="Pending images",
            listbox_height=15,
            dblclickcommand=self.dopendingtocurrent,
            )
        self.pendinglistbox.pack(side='top', padx=5, pady=5, fill='x')

        self.currentimageentry = Pmw.EntryField(fileframe,
            labelpos='n',
            label_text="Current image",
            )
        self.currentimageentry.pack(side='top', padx=5, pady=5, fill='x')

        self.finishedlistbox = Pmw.ScrolledListBox(fileframe,
            labelpos='n',
            label_text="Finished images",
            listbox_height=15,
            dblclickcommand=self.dofinishedtocurrent,
            )
        self.finishedlistbox.pack(side='top', padx=5, pady=5, fill='x')

        # minimap at bottom
        Tkinter.Label(fileframe, text="Navigation").pack(side='top', 
            padx=5, fill='x')
        self.minimap = MiniMap(fileframe, 175)
        self.minimap.pack(side='top', padx=5, pady=5, fill='x')

        self.minimap.addlistener(self.scrollto)


        # ----- main 2d view area
        viewframe = Tkinter.Frame(mainframe)
        viewframe.pack(side='left', expand='yes', fill='both', padx=5, pady=5)

        self.imageframe = Pmw.ScrolledFrame(viewframe,
            frame_cursor="crosshair",
            # tried and didn't like: tcross & cross (lot like crosshair), dotbox
            #   and target (way too thick), dot (opaque, too big)
            )
        self.imageframe.pack(side='top', expand='yes', fill='both')

        self._imagelabel = Tkinter.Label(self.imageframe.interior())
        self._imagelabel.pack(side='top', anchor='center')

        # interactions:
        self._imagelabel.bind("<Button-1>", self.doclickimage)
        self._imagelabel.bind("<Motion>", self.domousemotion)




        # ----- tools
        toolframe = Tkinter.Frame(mainframe)
        toolframe.pack(side='left', fill='both')


        # metadata area:
        metadatagroup = Pmw.Group(toolframe, 
            tag_text = "Metadata",
            )
        metadatagroup.pack(side='top', fill='x', padx=4, pady=2)

        self.tiltangleentry = Pmw.EntryField(metadatagroup.interior(),
            labelpos='w',
            label_text='tilt angle (deg):',
            validate={'validator': 'real'},
            modifiedcommand=self.dochangetiltangle,
            )
        self.tiltangleentry.pack(side='top', fill='x', padx=2, pady=2)

        # tools
        toolgroup = Pmw.Group(toolframe,
            tag_text="Tools",
            )
        toolgroup.pack(side='top', fill='x', padx=4, pady=2)

        self.toolradioselect = Pmw.RadioSelect(toolgroup.interior(),
            buttontype='radiobutton',
            orient='vertical',
            )
        self.toolradioselect.pack(side='top', fill='x', pady=2)
        for tool in toollist:
            self.toolradioselect.add(tool)
        self.toolradioselect.invoke(defaulttool)


        # points group
        pointgroup = Pmw.Group(toolframe, 
            tag_text = "Points",
            )
        pointgroup.pack(side='top', fill='x', padx=4, pady=2)




        # centering
        self.centeringlistbox = Pmw.ScrolledListBox(pointgroup.interior(),
            labelpos='n',
            label_text="Centering reflections",
            items=[],
            listbox_height = 6,
            )
        self.centeringlistbox.pack(side='top', fill='x', padx=2, pady=2)
        bbox = Pmw.ButtonBox(pointgroup.interior())
        bbox.pack(side='top', fill='x')
        bbox.add("Remove", command=self.doremovecenter)
        bbox.add("Clear", command=self.doclearcenters)


        # reference
        self.referencelistbox = Pmw.ScrolledListBox(pointgroup.interior(),
            labelpos='n',
            label_text="Reference points",
            items=[],
            listbox_height = 2,
            )
        self.referencelistbox.pack(side='top', fill='x', padx=2, pady=2)
        bbox = Pmw.ButtonBox(pointgroup.interior())
        bbox.pack(side='top', fill='x')
        bbox.add("Remove", command=self.doremoveref)
        bbox.add("Clear", command=self.doclearrefs)

        
        # beam stop                
        self.beamstopentry = Pmw.EntryField(pointgroup.interior(),
            labelpos='n',
            label_text="Beam stop center",
            )
        self.beamstopentry.pack(side='top', fill='x', padx=2, pady=2)

        # spots
        self.spotlistbox = Pmw.ScrolledListBox(pointgroup.interior(),
            labelpos='n',
            label_text="Spots",
            items=[],
            listbox_height = 10,
            )
        self.spotlistbox.pack(side='top', fill='x', padx=2, pady=2)

        spotframe = Tkinter.Frame(pointgroup.interior())
        spotframe.pack(side='top', fill='x')
        Tkinter.Checkbutton(spotframe, text="Refine", 
            variable=self.refinespotsvar).pack(side='left')
        self.shapemenu = Pmw.OptionMenu(spotframe,
            items=["circle", "square"],
            menubutton_width=6,
            initialitem="circle",
            )
        self.shapemenu.pack(side='left')
        self.shapesizeentry = Pmw.EntryField(spotframe,
            labelpos='e',
            label_text="pix",
            value=10,
            validate={'validator': 'numeric',
                'min': 1,
                'minstrict': 1,
                },
            entry_width=3,
            )
        self.shapesizeentry.pack(side='top', fill='x')

        bbox = Pmw.ButtonBox(pointgroup.interior())
        bbox.pack(side='top', fill='x')
        bbox.add("Remove", command=self.doremovespot)
        bbox.add("Clear", command=self.doclearspots)




        # ----- final preparations

        # clear initial data and UI
        self.clearimageanddata()
        self.updateall()

        # go!
        self.root.mainloop()
        
        # end __init__()
    
    # ......................... centerofmass() .........................
    def centerofmass(self, data):
        """
        given 2d numpy data array, find the center of mass
        
        input: 2d scalar array
        output: x, y center of mass (floating point)
        """
        
        ny, nx = data.shape

        # given some intensity field I(x, y), you want to do:
        # xcom = sum(xi * I(xi, yi)) / sum(I(xi, yi)), with the sum
        #   running over the coordinates xi and yi 
        # (likewise for ycom)

        xcoord = numpy.arange(nx, dtype=numpy.float)
        ycoord = numpy.arange(ny, dtype=numpy.float)
        xramp, yramp = numpy.meshgrid(xcoord, ycoord)

        datasum = data.sum()
        xc = (data * xramp).sum() / datasum
        yc = (data * yramp).sum() / datasum

        return xc, yc
        
        # end centerofmass()
        
    # ......................... clearimageanddata() .........................
    def clearimageanddata(self):
        """
        clear out the data; expected to be followed by a call to 
        updateall() to reset the UI as well
        """

        # reset image & data
        self.imagefilepath = ""
        self.image = blankimage
        self.imagedata = None
        self.datamax = -1

        # clear various points lists
        self._centers = []
        self._beamcenter = []
        self._references = []
        self._beamstopcenter = []
        self._spots = []
        self._tiltangle = 0.0

        # end clearimageanddata()
    
    # ......................... despeckle() .........................
    def despeckle(self, data):
        """
        3x3 median filter on the input data
        
        input: 2d numpy array
        output: despeckled 2d numpy array, same size, original unchanged
        """
        
        ny, nx = data.shape
        big = numpy.zeros((9, ny, nx), dtype=data.dtype)
        big[0] = data
        # y axis
        big[1] = numpy.roll(data, 1, axis=0)
        big[2] = numpy.roll(data, -1, axis=0)
        # x axis
        big[3] = numpy.roll(data, 1, axis=1)
        big[4] = numpy.roll(data, -1, axis=1)
        # up/down in y, then x
        big[5] = numpy.roll(big[1], 1, axis=1)
        big[6] = numpy.roll(big[1], -1, axis=1)
        big[7] = numpy.roll(big[2], 1, axis=1)
        big[8] = numpy.roll(big[2], -1, axis=1)

        return numpy.median(big, axis=0)

        # end despeckle()
        
    # ......................... dochangetiltangle() .........................
    def dochangetiltangle(self):
        """
        called when user changes value in tilt angle entry
        """
        
        self._tiltangle = self.tiltangleentry.getvalue()
        
        # end dochangetiltangle()
        
    # ......................... doclearbeamstop() .........................
    def doclearbeamstop(self):
        """
        clear beam stop entry
        """
        
        self._beamstopcenter = []
        self.updateall()
        # self.beamstopentry.clear()        
        
        # end doclearbeamstop()
                
    # ......................... doclearcenters() .........................
    def doclearcenters(self):
        """
        clear centers list
        """
        
        self._centers = []
        self._beamcenter = []
        self.updateall()
        
        # end doclearcenters()
                        
    # ......................... doclearrefs() .........................
    def doclearrefs(self):
        """
        clear reference spots
        """
        
        self._references = []
        self.updateall()
        
        # end doclearrefs()

    # ......................... doclearspots() .........................
    def doclearspots(self):
        """
        clear the spot list
        """

        self._spots = []
        self.updateall()
        
        # end doclearspots()
        
    # ......................... doclickimage() .........................
    def doclickimage(self, event):
        """
        called on left-click in image label
        """
        
        if not self.imagefilepath:
            # no image, do nothing
            return

        newspot = event.x, event.y
    
        tool = self.toolradioselect.getvalue()
        if tool == "centering reflections":
            self._centers.append(newspot)
        elif tool == "reference points":
            self._references.append(newspot)
        elif tool == "beam stop center":
            self._beamstopcenter = [newspot]
        elif tool == "spots":
            if self.refinespotsvar.get():
                newspot = self.recenter((event.x, event.y))
                if newspot is None:
                    tkMB.showerror(title="Error",
                        message="Couldn't refine spot!  Are you too close to the edge?")
                    return

                if self.options.debug:
                    print "click at %s; com at %s" % ((event.x, event.y), newspot)
            self._spots.append(newspot)

        self.updateall()
        
        # end doclickimage()
        
    # ......................... docurrenttofinished() .........................
    def docurrenttofinished(self):
        """
        move current image to finished list
        """
        
        filename = self.currentimageentry.getvalue()
        if not filename:
            return
        
        # save info
        self._data["images"][filename] = self.getcurrentdata()
        self.save()

        self.clearimageanddata()
        self.updateall()
        
        self.currentimageentry.setvalue("")

        items = list(self.finishedlistbox.get())
        items.append(filename)
        self.finishedlistbox.setlist(sorted(items))

        # end docurrenttofinished()
        
    # ......................... dofileopendirectory() .........................
    def dofileopendirectory(self):
        """
        open a directory's images
        """
        
        # if already dir open, save
        if len(self._data) > 0:
            self.save()


        # prompt for dir
        self._imagedirectory = tkFD.askdirectory(title="Open directory...")
        if not self._imagedirectory:
            return


        # read dir; if save file exists, load into data structs;
        savefilepath = os.path.join(self._imagedirectory, savefilename)
        if os.path.exists(savefilepath):
            tempdata = json.load(open(savefilepath))
            if tempdata["metadata"]["description"] != savefiledescription:
                tkMB.showerror(title="No!",
                    message="Despite the name, %s doesn't appear to be a cataspot data file!" %
                    savefilepath)
                return

            if tempdata["metadata"]["file version"] > savefileversion:
                tkMB.showerror(title="No!",
                    message="=%s is from a newer version of cataspot!")
                return

            self._data = tempdata["data"]
            self.launderdata()

            # populate pending/finished lists; this is a bit tricky; we
            #   only populate files that are present now; if anything's
            #   been added since last save, it's shown now; if anything's
            #   missing, it's kept in the saved data, but it's not shown
            # we might want to adjust this later
            imagefilelist = self.loadnewdirectory(self._imagedirectory)
            finishedfilelist = [filename for filename in imagefilelist 
                if filename in self._data["images"]]
            pendingfilelist = [filename for filename in imagefilelist
                if filename not in finishedfilelist]

        else:
            # no existing data, start anew
            self._data = {}
            self._data["directory"] = self._imagedirectory
            self._data["images"] = {}

            pendingfilelist = self.loadnewdirectory(self._imagedirectory)
            finishedfilelist = []

        self.pendinglistbox.setlist(sorted(pendingfilelist))
        self.finishedlistbox.setlist(sorted(finishedfilelist))
        
        self.folderentry.setvalue(os.path.basename(self._imagedirectory))
        self.clearimageanddata()
        self.updateall()

        # end dofileopendirectory()
        
    # ......................... dofilequit() .........................
    def dofilequit(self):
        """
        File > Quit
        """
        
        # autosave before quit
        self.save()

        self.root.destroy()
        
        # end dofilequit()
        
    # ......................... dofinishedtocurrent() .........................
    def dofinishedtocurrent(self):
        """
        move clicked-on item from finished list to current
        """
        
        # get item
        filename = self.finishedlistbox.getcurselection()[0]

        # move any current to finished
        if self.currentimageentry.getvalue():
            self.docurrenttofinished()

        # move finished to current
        self.setcurrentimage(filename)

        # remove from finished list:
        items = list(self.finishedlistbox.get())
        items.remove(filename)
        self.finishedlistbox.setlist(sorted(items))
        
        # end dofinishedtocurrent()
        
    # ......................... dohelpabout() .........................
    def dohelpabout(self):
        """
        Help > About
        """
        
        tkMB.showinfo(title="About",
            message="Cataspot v%s\nby Donald J. Olbris" % __version__)
        
        # end dohelpabout()
        
    # ......................... domousemotion() .........................
    def domousemotion(self, event):
        """
        called when mouse moves
        """
        
        x, y = event.x, event.y
        if self.imagedata is None:
            value = 0
        else:
            value = self.imagedata[y, x]
        self.coordinatelabel.configure(text=coordinatetemplate % (x, y, value))
        
        # end domousemotion()
        
    # ......................... dopendingtocurrent() .........................
    def dopendingtocurrent(self):
        """
        move clicked-on pending item to current
        """
        
        # get item
        filename = self.pendinglistbox.getcurselection()[0]

        # move any current to finished
        if self.currentimageentry.getvalue():
            self.docurrenttofinished()

        # move pending to current
        self.setcurrentimage(filename)

        # remove from pending list:
        items = list(self.pendinglistbox.get())
        items.remove(filename)
        self.pendinglistbox.setlist(sorted(items))
        
        # end dopendingtocurrent()
        
    # ......................... doremovecenter() .........................
    def doremovecenter(self):
        """
        remove a center
        """
        
        ans = self.centeringlistbox.getvalue()
        if ans:
            item = ans[0]
            self._centers.remove(item)
            self.updateall()
        
        # end doremovecenter()
        
    # ......................... doremoveref() .........................
    def doremoveref(self):
        """
        remove the selected reference
        """
        
        ans = self.referencelistbox.getvalue()
        if ans:
            item = ans[0]
            self._references.remove(item)
            self.updateall()
        
        # end doremoveref()
        
    # ......................... doremovespot() .........................
    def doremovespot(self):
        """
        remove the selected spot
        """
        
        ans = self.spotlistbox.getvalue()
        if ans:
            item = ans[0]
            self._spots.remove(item)
            self.updateall()
        
        # end doremovespot()
        
    # ......................... drawmarkers() .........................
    def drawmarkers(self, im, pointlist, marker):
        """
        draws a markers on the image
        
        input: PIL image; list of (x, y) locations; Marker instance
        output: none (draws on image)
        """
        
        if pointlist:
            draw = ImageDraw.Draw(im)
            for x, y in pointlist:
                if marker.kind == "x":
                    draw.line([x - marker.size, y - marker.size, 
                        x + marker.size, y + marker.size], 
                        fill=marker.color)
                    draw.line([x - marker.size, y + marker.size, 
                        x + marker.size, y - marker.size], 
                        fill=marker.color)
                elif marker.kind == "+":
                    draw.line([x - marker.size, y, x + marker.size, y],
                        fill=marker.color)
                    draw.line([x, y + marker.size, x, y - marker.size],
                        fill=marker.color)
                elif marker.kind == "o":
                    # draw.arc can't handle floats for some reason (line can);
                    #   so shift x, y, to nearest int for this marker
                    x = int(round(x))
                    y = int(round(y))
                    draw.arc([x - marker.size, y - marker.size, 
                        x + marker.size, y + marker.size], 
                        0, 360, fill=marker.color)
                elif marker.kind == "(+)":
                    # (+) = circle with crosshairs
                    # draw.arc can't handle floats for some reason (line can);
                    #   so shift x, y, to nearest int for this marker
                    x = int(round(x))
                    y = int(round(y))
                    draw.arc([x - marker.size, y - marker.size, 
                        x + marker.size, y + marker.size], 
                        0, 360, fill=marker.color)
                    draw.line([x - marker.size, y, x + marker.size, y],
                        fill=marker.color)
                    draw.line([x, y + marker.size, x, y - marker.size],
                        fill=marker.color)
                else:
                    raise ValueError("unknown marker kind %s" % marker.kind)
                
        # end drawmarkers()
        
    # ......................... getcurrentdata() .........................
    def getcurrentdata(self):
        """
        returns the internal data as a dictionary
        """

        datadict = {}
        datadict["centers"] = self._centers
        datadict["beamcenter"] = self._beamcenter
        datadict["references"] = self._references
        datadict["beamstopcenter"] = self._beamstopcenter
        datadict["spots"] = self._spots
        datadict["tiltangle"] = self._tiltangle

        return datadict

        # end getcurrentdata()
           
    # ......................... getmetadata() .........................
    def getmetadata(self):
        """
        get a basic metadata dictionary for json file saves
        """
        
        d = {}
        d["file version"] = savefileversion
        d["date"] = time.asctime()
        d["username"] = getpass.getuser()
        return d
                
        # end getmetadata()

    # ......................... launderdata() .........................
    def launderdata(self):
        """
        our points turn from tuples to lists when saved/loaded as json,
        so turn them back; otherwise they will behave badly in our
        list boxes, which inexplicably maintain tuple identity, but
        turn lists into strings that eval to lists, when you put 
        them into then take out of the list
        """
        
        # I've got to find a more elegant way to handle all these 
        #   lists...they all behave much the same, so I shouldn't
        #   be copying/pasting all this code...
        for image in self._data["images"]:
            for key in [
                "centers",
                "references",
                "beamstopcenter",
                "spots",
                ]:
                self._data["images"][image][key] = [tuple(item) for item in self._data["images"][image][key]]

        # end launderdata()
        
    # ......................... loadnewdirectory() .........................
    def loadnewdirectory(self, directorypath):
        """
        load a new directory; create new empty data structures
        
        input: path to directory
        output: list of filenames (just names, no paths)
        """
        
        # filter out images:
        results = []
        for fn in os.listdir(directorypath):
            base, extension = os.path.splitext(fn)
            if extension.lower() in [".tif"]:
                results.append(fn)
        return results

        # end loadnewdirectory()
        
    # ......................... openimage() .........................
    def openimage(self, filename):
        """
        open an image
        """
        
        with WatchCursorDisplayed(self.imageframe.component("frame")):
            try:
                image = Image.open(filename)
                filepath = os.path.abspath(filename)
            except:
                tkMB.showerror(title="Error!", 
                    message="There was an error reading %s!" % filename)
                return

            # worked this out in smvtools; expect it applies now, too
            if image.mode not in ["I;16", "I;16B"]:
                tkMB.showerror(title="Error!", 
                    message="expected 16-bit integers, got %s" % image.mode)
                return

            # clear data and start working
            self.clearimageanddata()

            self.image = image
            self.imagefilepath = filepath

            # set minimap; comes before despeckle, etc.
            self.minimap.setimage(self.image)

            self.imagedata = numpy.array(self.image.getdata(), dtype=numpy.uint16)
            nx, ny = self.image.size
            self.imagedata = self.imagedata.reshape(ny, nx)

            # experimental despeckle = 3x3 median filter: I don't want to
            #   add a scipy dependency, so brute force it; take image and
            #   shift by one pixel up and down in x and/or y eight times,
            #   assemble into a 3D array, and take the median on the z-axis;
            #   ignore edge effects where the shift wraps around, as 
            #   nothing we care about happens at the edges
            if self.despecklevar.get():
                start = time.time()
                self.imagedata = self.despeckle(self.imagedata)
                end = time.time()
                if self.options.debug:
                    print "despeckle: %ss" % (end - start)


            # figure out contrast parameters/colormap
            self.datamax = self.imagedata.max()


            # figure out scaling (not yet)


            # set title bar
            self.root.title("%s - %s" % (defaulttitle, os.path.basename(self.imagefilepath)))

            # update 
            self.updateall()

            # center the scrolled frame; update idle tasks is required so
            #   the scroll bars update after the image load, as their 
            #   position is used during the scrollto() operation
            self.imageframe.update_idletasks()
            ny, nx = self.imagedata.shape
            self.scrollto((nx // 2, ny // 2))

        # end openimage()
        
    # ......................... recenter() .........................
    def recenter(self, point):
        """
        recenter the input point
        
        input: (x, y)
        output: new (x, y)
        """
        
        xc, yc = point

        shape = self.shapemenu.getvalue()
        halfwidth = int(self.shapesizeentry.getvalue())
        width = 2 * halfwidth + 1
        ny, nx = self.imagedata.shape

        # extract data (copy)
        x1, y1 = xc - halfwidth, yc - halfwidth
        x2, y2 = xc + halfwidth, yc + halfwidth

        # can't be too close to edge:
        if x1 < 0 or y1 < 0 or x2 >= nx or y2 >= ny:
            return None

        # find COM in extract
        region = self.imagedata[y1: y2 + 1, x1: x2 + 1]

        if shape == "square":
            # no mask
            mask = 1
        elif shape == "circle":
            xcoord = numpy.arange(-halfwidth, halfwidth + 1, dtype=numpy.float)
            ycoord = numpy.arange(-halfwidth, halfwidth + 1, dtype=numpy.float)
            xramp, yramp = numpy.meshgrid(xcoord, ycoord)
            rsquared = xramp * xramp + yramp * yramp
            mask = numpy.where(rsquared <= halfwidth, 1, 0)

        localcenter = self.centerofmass(mask * region) 

        # re-offset point to image coords
        lx, ly = localcenter

        # NOTE: not clear if we want fp or int yet! 
        #   use int for now (easier to test, display in list);
        #   switch later if needed
        # lx = int(round(lx))
        # ly = int(round(ly))

        return lx + x1, ly + y1
        
        # end recenter()
        
    # ......................... save() .........................
    def save(self):
        """
        save the data 
        """
        
        if self._data:
            # throw current image into data, too, for save purposes:
            filename = self.currentimageentry.getvalue()
            if filename:
                self._data["images"][filename] = self.getcurrentdata()

            # I'm following what I do in Fly EM for JSON files re: 
            #   data and metadata
            d = {}
            d["metadata"] = self.getmetadata()
            d["metadata"]["description"] = savefiledescription

            d["data"] = self._data

            # somewhat amusingly, the "loading" image ends up in the data...
            if "loading" in d["data"]["images"]:
                del d["data"]["images"]["loading"]

            with open(os.path.join(self._imagedirectory, savefilename), 'wt') as f:
                json.dump(d, f, indent=2)

        # end save()

    # ......................... scrollto() .........................
    def scrollto(self, loc):
        """
        scroll the image so the given location is at the center of
        the window, or your best try
        
        input: (x, y) location
        output: none (scrolls frame)
        """
        
        if not self.imagefilepath:
            return

        # x/yview with "moveto" want the fraction of the window
        #   that should be offscreen to the left, which is about
        #   a ridiculous piece of input you could imagine

        x, y = loc
        ny, nx = self.imagedata.shape

        x1, x2 = self.imageframe.xview()
        size = x2 - x1
        if x < 0:
            xtarget = 0
        elif x >= nx - 1:
            xtarget = 1 - size
        else:
            xtarget = float(x) / nx - size / 2
        if xtarget < 0:
            xtarget = 0
        if xtarget > 1:
            xtarget = 1

        y1, y2 = self.imageframe.yview()
        size = y2 - y1
        if y < 0:
            ytarget = 0
        elif y >= ny - 1:
            ytarget = 1 - size
        else:
            ytarget = float(y) / ny - size / 2
        if ytarget < 0:
            ytarget = 0
        if ytarget > 1:
            ytarget = 1

        self.imageframe.xview('moveto', xtarget)
        self.imageframe.yview('moveto', ytarget)
        
        # end scrollto()
        
    # ......................... setcurrentdata() .........................
    def setcurrentdata(self, datadict):
        """
        set the internal data to the contents of the dictionary
        """

        self._centers = datadict["centers"]
        self._beamcenter = datadict["beamcenter"]
        self._references = datadict["references"]
        self._beamstopcenter = datadict["beamstopcenter"]
        self._spots = datadict["spots"]
        self._tiltangle = datadict["tiltangle"]

        # end setcurrentdata()
                    
    # ......................... setcurrentimage() .........................
    def setcurrentimage(self, filename):
        """
        load image and its data, if it's got any
        """
            
        # load image
        self.currentimageentry.setvalue("loading")
        self.openimage(os.path.join(self._imagedirectory, filename))

        # set label; remove from pendinglist
        self.currentimageentry.setvalue(filename)

        # check existing data; if we've got any, load it into the UI
        if filename in self._data["images"]:
            self.setcurrentdata(self._data["images"][filename])
        self.updateall()
        
        # end setcurrentimage()
        
    # ......................... _setimagelabel() .........................
    def _setimagelabel(self, image):
        """
        put the image into the image label
        
        input: Tk image
        """
        
        self._tkimage = ImageTk.PhotoImage(image)
        self._imagelabel.configure(image=self._tkimage)
        
        # end _setimagelabel()
        
    # ......................... setupmenus() .........................
    def setupmenus(self):
        """
        set up main menu
        """
        
        self.menubar = Pmw.MainMenuBar(self.root)
        self.root.configure(menu=self.menubar)
        
        # File
        self.menubar.addmenu('File', "")
        self.menubar.addmenuitem('File', 'command',
            command=self.dofileopendirectory, label="Open directory...")
        self.menubar.addmenuitem('File', 'command',
            command=self.dofilequit, label="Quit")
        

        # Options
        self.menubar.addmenu('Options', "")
        self.despecklevar = Tkinter.IntVar()
        self.despecklevar.set(1)
        self.menubar.addmenuitem("Options", "checkbutton",
            variable=self.despecklevar,
            label="Despeckle images")
        self.refinespotsvar = Tkinter.IntVar()
        self.refinespotsvar.set(1)
        self.menubar.addmenuitem("Options", "checkbutton",
            variable=self.refinespotsvar,
            label="Refine spot locations")

        # Help
        self.menubar.addmenu('Help', "")
        self.menubar.addmenuitem('Help', 'command',
            command=self.dohelpabout, label="About...")
        
        # end setupmenus()

    # ......................... updateall() .........................
    def updateall(self):
        """
        update all GUI stuff
        """
        
        self.updatebeamstopcenter()
        self.updatecenterlist()
        self.updatereflist()
        self.updatespotlist()
        self.updatetiltangle()
        self.updateimage()
        
        # end updateall()
        
    # ......................... updateimage() .........................
    def updateimage(self):
        """
        update the image
        """
        
        # grab the image
        # self._displayimage = self.image.convert('RGB')

        # contrast: scale so max is about 10% of data max, and scale down to 8-bit
        if self.imagedata is not None:
            scalefactor = self.datamax * maxcontrast
            showdata = (255. / scalefactor) * self.imagedata.clip(max=scalefactor)
            self._displayimage = Image.fromarray(showdata.astype(numpy.uint8))
            self._displayimage = self._displayimage.convert('RGB')
        else:
            # no data, so grab the default blank image
            self._displayimage = blankimage

        # draw circles on spots:
        self.drawmarkers(self._displayimage, self._centers, toolmarkerdict["centering reflections"])
        self.drawmarkers(self._displayimage, self._beamcenter, toolmarkerdict["beam center"])
        self.drawmarkers(self._displayimage, self._references, toolmarkerdict["reference points"])
        self.drawmarkers(self._displayimage, self._beamstopcenter, toolmarkerdict["beam stop center"])
        self.drawmarkers(self._displayimage, self._spots, toolmarkerdict["spots"])


        # scale
        # (skip for now)




        # put the Tk into the label
        self._setimagelabel(self._displayimage)

        
        # end updateimage()
    
    # ......................... updatebeamstopcenter() .........................
    def updatebeamstopcenter(self):
        """
        update beam stop center entry from internal data
        """
    
        if self._beamstopcenter:        
            self.beamstopentry.setvalue(self._beamstopcenter[0])
        else:
            self.beamstopentry.clear()
        
        # end updatebeamstopcenter()
        
    # ......................... updatecenterlist() .........................
    def updatecenterlist(self):
        """
        update list of centers from internal data; calculate beam center
        """
        
        self.centeringlistbox.setlist(sorted(self._centers))

        if self._centers:
            xsum = sum(x for x, y in self._centers)
            ysum = sum(y for x, y in self._centers)
            self._beamcenter = [(float(xsum) / len(self._centers), float(ysum) / len(self._centers))]
        else:
            self._beamcenter = []
        
        # end updatecenterlist()
    
    # ......................... updatereflist() .........................
    def updatereflist(self):
        """
        update list of reference from internal data
        """
        
        self.referencelistbox.setlist(sorted(self._references))
        
        # end updatereflist()
                        
    # ......................... updatespotlist() .........................
    def updatespotlist(self):
        """
        update the list of spots from the internal data
        """
        
        self.spotlistbox.setlist(sorted(self._spots))
        
        # end updatespotlist()

    # ......................... updatetiltangle() .........................
    def updatetiltangle(self):
        """
        update tilt angle from internal data
        """
        
        self.tiltangleentry.setvalue(self._tiltangle)
        
        # end updatetiltangle()
        
    # end class CataspotApp

# ------------------------- class MiniMap -------------------------
class MiniMap(Tkinter.Frame):
    """
    this class is a minimap for Tkinter
    """
    # ......................... __init__ .........................
    def __init__(self, parent, width, **kwargs):
        """
        
        input:  parent window
                width of minimap (height determined by aspect ratio of
                    input image)
        """
        
        # superclass init:
        Tkinter.Frame.__init__(self, parent, **kwargs)        
        self.width = width


        # the GUI itself is fairly simple: a label holding an image
        self._imagelabel = Tkinter.Label(self)
        self._imagelabel.pack()

        # initial image is blank, exact size, and square:
        im = Image.new("RGB", (width, width))
        self.setimage(im)

        # interaction
        self._imagelabel.bind("<Button-1>", self.doclickimage)

        # callbacks
        self._listeners = []

        # end __init__()
    
    # ......................... addlistener() .........................
    def addlistener(self, listener):
        """
        add listener for clicks
        
        input: callable taking (x, y) tuple in image coordinates
        output: none
        """
        
        self._listeners.append(listener)
        
        # end addlistener()
        
    # ......................... doclickimage() .........................
    def doclickimage(self, event):
        """
        called on click in image; triggers callbacks
        """
        
        imageloc = int(event.x / self._scalefactor), int(event.y / self._scalefactor)

        for listener in self._listeners:
            listener(imageloc)
    
        # end doclickimage()

    # ......................... removelistener() .........................
    def removelistener(self, listener):
        """
        removes a previously added listener; silent if it doesn't exist!
        """
        
        if listener in self._listeners:
            self._listeners.remove(listener)
        
        # end removelistener()
        
    # ......................... setimage() .........................
    def setimage(self, im):
        """
        set the image for the minimap
        """
        
        # first, scale it down (never up)
        self._originalimage = im
        nx, ny = im.size
        self._originalsize = nx, ny

        self._scalefactor = float(self.width) / nx
        if nx > self.width:
            height = int(self._scalefactor * ny)
            im = im.resize((self.width, height))
        else:
            height = ny
        self._scaledsize = self.width, height
        self._image = im
        self._tkimage = ImageTk.PhotoImage(self._image)
        self._imagelabel.configure(image=self._tkimage)
        
        # end setimage()
        
    
    
    # end class MiniMap

# ------------------------- class WatchCursorDisplayed -------------------------
class WatchCursorDisplayed(object):
    """
    context manager for showing watch cursor; intended to be used like this:

    with WatchCursorDisplayed(parentwindow):
        # do stuff

    """
    # ......................... __init__ .........................
    def __init__(self, parent):
        """
        
        """
        
        # store values for now

        self.parent = parent
        
        # end __init__()
    
    # ......................... __enter__() .........................
    def __enter__(self):
        """
        called when context manager is entered
        """

        self.oldcursor = self.parent.cget("cursor")
        self.parent.configure(cursor="watch")
        self.parent.update()
        return self
        
        # end __enter__()
                
    # ......................... __exit__() .........................
    def __exit__(self, type, value, tb):
        """
        called when context manager exits
        
        input: type, value, traceback as returned by sys.exc_info()
        output: none
        """
        
        self.parent.configure(cursor=self.oldcursor)
        self.parent.update()
        
        # end __exit__()
        
    # end class WatchCursorDisplayed

# ------------------------- main() -------------------------
def main():
    """
    run the program!
    """
    
    usage = "usage: %prog [options]"
    parser = optparse.OptionParser(usage=usage, version="%%prog %s" % __version__)
    parser.add_option("-d", "--debug", action="store_true", dest="debug", default=False,
        help="enable debug mode")
        
    (options, args) = parser.parse_args()

    CataspotApp(options)
    
    # end main()

if __name__ == '__main__':
    main()