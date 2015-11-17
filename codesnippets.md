---
layout: page
title: Code Snippets
permalink: /Code Snippets/
---

Here, I list some code snippets from previous projects I've worked on. 

The following is from some of the work I did at the MIIL. The function takes in 
data from a PET (Positron Emission Tomography) scanner and converts it into a 
visual representation of the source of radiation. This function was used to
diagnose problems with the system and visualize the data collected.

{% highlight python %}
import pylab as pl
import numpy as np
import matplotlib.pyplot as plt
from math import sqrt
import time

def simple_bp(f_lst, disp=True, mslope=None, first_cathode=False):
    '''
    --Description--
    Outputs simple unfiltered backprojection
    with accompanying system matrix
    --Input Details--
    f_lst: *.lst file to reconstruct
    disp: whether or not to display
    mslope: minimum permitted slope in reconstruction
    first_cathode: whether or not to include only first cathode events
    '''
    start = time.time()
    ##### INITIALIZATION #####
    x1y1 = []
    x2y2 = []#coordinates of coincidence events
    with open(f_lst, 'r') as f:
        for line in f:
            l = map(float, line.split())
            if first_cathode:
                c1 = ' '.join(line.split()[:4])
                c2 = ' '.join(line.split()[9:13])
                if electrode(c1).split()[1] != '1' or electrode(c2).split()[1] != '1':
                    continue
                    #ignofre if not first cathode
            x1y1 += [(l[6], l[7])] #(x, y) coordinate pair
            x2y2 += [(l[15], l[16])]
    ##### we have finished reading in our coordinates #####
    for i, xy in enumerate(x1y1): #make all coordinates postive for matrix indices
        x1y1[i] = [xy[0] - X_MIN, xy[1] - Y_MIN] #the min will be negative 
    for i, xy in enumerate(x2y2):
        x2y2[i] = [xy[0] - X_MIN, xy[1] - Y_MIN]
    ##### Initializing the system matrix #####
    system = np.zeros((140 * SCALE, 100 * SCALE ))
    ##### FINISHED INITIALIZING #####
     
    #The system matrix will be a matrix of integers describing the magnitude of
    #LOR density in that particular pixel. This essentially shows the density of
    #radioisotopes in the sample tested.
    # IMPORTANT: (Row, Column) => (y, x)
    LORs = zip(x1y1, x2y2) #zipping together endpoints of LOR
    print len(LORs)
    if not mslope:
        for LOR in LORs:
            backproject(LOR, system) #do backprojections on system matrix
    else:
        for LOR in LORs:
            backproject(LOR, system, min_slope=mslope)
    ##### Displaying System Matrix #####
    #cmap='hot'#
    print 'Done, that took %f seconds'%(time.time()-start)
    if disp:
        if mslope:
            pl.title('Minimum Slope Permitted: %d'%(mslope))
        pl.imshow(system)
        pl.colorbar()
        pl.show()
         
        return system
    else:
        return system
{% endhighlight %}

![Alt text](/images/dls_recon.png "Sample Output of the Above Function")

Here is a sample of some of my GUI code from the same project. This class allows matplotlib graphs
to be displayed in the GUI. The GUI itself provides testing and debugging procedures for the PET system. 

{% highlight python %}

class MPLPanel(wx.Panel):
    def __init__(self, parent):
        wx.Panel.__init__(self, parent, -1, size=(500, 600))
        self.figure = Figure()
        gs = gridspec.GridSpec(2, 1, height_ratios=[3, 1])
        self.axes_dist = self.figure.add_subplot(gs[0])
        self.axes_flood = self.figure.add_subplot(gs[1])
        self.DataAcquired = False
        self.FHAcquired = False
        self.canvas = FigureCanvas(self, -1, self.figure)
        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(self.canvas, 1, wx.EXPAND|wx.ALL)
        self.SetSizerAndFit(sizer)
        #self.canvas.Centre()
    def graph(self, NBRC, f_unpacked=None, f_calib=None, numLines=2 * 10 ** 7, kev=False):
        #self.figure.clf(keep_observers=True)
        if not self.DataAcquired:
            print 'ACQUIRING DATA'
            self.data = energyhist.energydata(f_unpacked, f_calib, numLines)
            self.DataAcquired = True
            print 'DATA ACQUIRED'
        if kev:
            self.axes_dist.hist(self.data[NBRC][0],bins=200, histtype='step', stacked=True, fill=True)
            self.axes_dist.set_xlabel('Energy (ADC)')
        else:
            self.axes_dist.hist(self.data[NBRC][1],bins=200, histtype='step', stacked=True, fill=True)
            self.axes_dist.set_xlabel('Energy (keV)')
        self.axes_dist.set_title('Energy Distribution')
        self.axes_dist.set_ylabel('Counts')
        self.canvas = FigureCanvas(self, -1, self.figure)
        print 'ENERGY DISTRIBUTION GRAPHED'
    def flood_hist(self, node, board, f_acp=None):
        if not self.FHAcquired:
            print 'ACQUIRING DATA'
            self.fhdata = activity(f_acp)
            self.FHAcquired = True
            print 'DATA ACQUIRED'
        self.axes_flood.imshow(self.fhdata[node, board])
        #self.axes_flood.set_title('Flood Histogram of Coincidence Pairs')
        self.axes_flood.set_xlabel('Anode Number')
        self.axes_flood.set_ylabel('Cathode Number')
        self.canvas = FigureCanvas(self, -1, self.figure)
        print 'FLOOD HISTOGRAM GENERATED'

{% endhighlight %}

![Alt text](/images/GUIPVCZT.png "Screenshot of the Full GUI")

Here's python code I wrote (just for fun) that outputs all permutations of a list (of distinct elements
-- we could handle lists that have certain identical elements with a simple conditional in the middle).

{% highlight python %}

def permutations(lst):
    perms = []
    if len(lst) <= 1:
        perms += [lst]
    else:
        for i in range(len(lst)):
            temp = lst[i]
            lst[i] = lst[0]
            lst[0] = temp
            perms += [[lst[0]] + p for p in permutations(lst[1:])]
    return perms

{% endhighlight %}

Here's python code I wrote that generates all subsets of a given list:

{% highlight python %}

def subsets(lst):
    if len(lst) == 0:
        return [lst]
    else:
        cdrsubs = subsets(lst[1:])
        return cdrsubs + [[lst[0]] + s for s in cdrsubs]

{% endhighlight %}
From my adventures in functional programming... This version of quicksort, written in Scheme (*shudders*
... so many parentheses ...) is inspired by the version featured
on the Haskell website when I was exploring that language. (Assume filter has previously been defined)

{% highlight scheme %}
(define (qsort lyst)
  (if (null? lyst)
      lyst
      (append (qsort (filter (lambda (x) (<= x (car lyst))) (cdr lyst)))
              (cons (car lyst)
                    (qsort (filter (lambda (x) (> x (car lyst))) (cdr lyst)))))))
{% endhighlight %}
