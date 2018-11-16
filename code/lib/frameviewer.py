# $Author: patricio $
# $Revision: 285 $
# $Date: 2010-06-18 17:59:25 -0400 (Fri, 18 Jun 2010) $
# $HeadURL: file:///home/esp01/svn/code/python/branches/patricio/photpipe/lib/frameviewer.py $
# $Id: frameviewer.py 285 2010-06-18 21:59:25Z patricio $

#! /usr/bin/env python

import os, sys
import pyfits as pf
import numpy  as np
import matplotlib.pyplot as plt
import pygame            as pg
from   pygame.locals import *
from   manageevent   import *

def frameviewer(event, zoom=True, expand=14):
  """
    Display in screen the masked data frames, opverplotting the position 
    from the centering routine, and aperture sizes (is exist).

    Parameters
    ----------
    event : An Event instance
            An event after poet_3center. If an event after poet_4photom is
            is supplied, it  will display the star and sky apertures as 
            well.
    zoom  : boolean
            If True, display a zoom area around the star. Necessary for 
            large arrays.

    Example
    -------
    >>> event = loadevent('tr001bs51_pht')
    >>> from frameviewer import *
    >>> frameviewer(event, zoom=True)

    Usage
    -----
    Arrow keys:
    Left / Right: move through different positions.
    Up / Down   : move through different frames in the same position.

    To go to a certain frame, type the number of a frame and then press 
    Enter (don't use the numeric keypad).

    Revisions
    ---------
    2010-07-12  Written by: Patricio Cubillos       pcubillos@fulbrightmail.org
    2010-07-13  patricio    Added documentation.
  """
  
  # look for data:
  try:
    dummy = event.data[0,0,0,0]
  except:
    curdir = os.getcwd()
    os.chdir(event.eventdir + '/' + event.centerdir)
    updateevent(event, event.eventname + "_ctr", ['data', 'mask'])
    os.chdir(curdir)

  # Look for phtometry apertures:
  try:
    photap = event.photap
    skyin  = event.skyin
    skyout = event.skyout
    aper = True
  except:
    aper = False

  maxnimpos, ny, nx, npos = np.shape(event.data)

  realny, realnx = ny, nx  # keep the values when we zoom in

  if zoom:   # zoom in, around the star.
    expand = 20
    ny, nx = 35, 35

  height, width = ny*expand, nx*expand
  fps = 10 # Frames per second

  realscreen = pg.display.set_mode((width, height)) # The real screen
  image      = pg.Surface((width, height))

  red       = 255,   0,   0
  blue      =   0,   0, 255
  green     =   0, 255,   0
  black     =   0,   0,   0
  white     = 255, 255, 255
  maskcolor = 255,  70,   0
  oob       =  25,  25, 112

  clock  = pg.time.Clock() #The clock for controlling FPS

  pg.init()
  font = pg.font.SysFont(None, 20)  # 20 = size

  def drawcross(surf, x, y, col, len=3.0, open=2.0):
    pg.draw.line(surf, col, (x-len-open, y), (x-open,    y),  1)
    pg.draw.line(surf, col, (x+open,     y), (x+len+open,y),  1)
    pg.draw.line(surf, col, (x, y-len-open), (x, y-open),     1)
    pg.draw.line(surf, col, (x, y+open),     (x, y+len+open), 1)

  def plot_text(font, texto, screen, position, color):
    the_text = font.render( texto.strip(), True, color )
    rect = the_text.get_rect()
    rect.midtop = position
    screen.blit(the_text, rect)


  # load first image
  def loadimage(index, posit):

    if zoom:
      im = np.zeros((ny, nx))
      ma = -np.ones((ny, nx))
      center = event.srcest[:,posit]
      uplim = center[0]+ny/2+1
      lolim = center[0]-ny/2
      rilim = center[1]+nx/2+1
      lelim = center[1]-nx/2

      for   i in np.arange(np.amax((0,lelim)), np.amin((realnx,rilim))):
        for j in np.arange(np.amax((0,lolim)), np.amin((realny,uplim))):
          im[j-lolim, i-lelim] = event.data[index, j, i, posit]
          ma[j-lolim, i-lelim] = event.mask[index, j, i, posit]

    else:
      im = event.data[index,:,:,posit]
      ma = event.mask[index,:,:,posit]

    # mask bad pixels
    mim = np.copy(im)
    mim[np.where(ma==0)] = 0

    # Normalize colors to range [0,255]
    imax, imin = np.amax(mim), np.amin(mim)
    mim = (mim - imin) * 255 / (imax - imin)

    # Draw the image
    for   i in np.arange(nx):
      for j in np.arange(ny):
        if ma[-j-1,i] == 1:  # if the pixel is good draw it
          image.fill( (mim[-j-1,i], mim[-j-1,i], mim[-j-1,i], 1.0), 
                      (i*expand, j*expand, (i+1)*expand, (j+1)*expand) )
        elif ma[-j-1,i] == -1: # out of bounds (for zoom == True)
          image.fill( oob, 
                      (i*expand, j*expand, (i+1)*expand, (j+1)*expand))
        else:           # Draw the mask
          image.fill( maskcolor, 
                      (i*expand, j*expand, (i+1)*expand, (j+1)*expand))
  
    # Draw the center
#    print(' ' + str(event.fp.x[index,posit])+' '+str(event.fp.y[index,posit]))
    y, x = ny-event.fp.y[posit,index]-0.5, event.fp.x[posit,index]+0.5
    if zoom:
      x =  event.fp.x[posit,index] - event.srcest[1,posit] + nx/2 + 0.5
      y = -event.fp.y[posit,index] + event.srcest[0,posit]   + ny/2 + 0.5    

    drawcross(image, x*expand, y*expand, red)

    if aper:
      pg.draw.circle(image, red,  (x*expand, y*expand), photap*expand, 1)
      pg.draw.circle(image, blue, (x*expand, y*expand), skyin *expand, 1)
      pg.draw.circle(image, blue, (x*expand, y*expand), skyout*expand, 1)


  # Captions
  caption = event.planetname + ' Image IM in position POS'
  pg.display.set_caption(caption) #Set the caption to the first one
  level = 1


  def game(event):

    id   = 0
    pos  = 0
    goto = ''

#    global id
#    global pos  #So we can change the caption vars
    loadimage(id,pos)

    while True:
      pg.display.set_caption( caption.replace('IM', str(id)).replace('POS', str(pos)) )
      
      for event in pg.event.get(): #The all important event loop
        if event.type == QUIT: 
          pg.quit()
          sys.exit() # press the X in the corner
  
        if event.type == KEYDOWN: #If user pressed a key
          if event.key == K_ESCAPE: 
            pg.quit()
            sys.exit() #Exit is pressed Esc
          if event.key == K_F1: #F1 key
            #Screen shot.
            for i in range(9999): #See if file mow 0000.bmp exists.
              path = os.path.join(os.getcwd(),'mow '+str(i).zfill(4)+'.bmp')
              if os.path.exists(path): #if not, see if mow 0001.bmp etc...
                continue
              else: #If not found, save it as that filename
                #If shift is held, save the screen
                if pg.key.get_pressed()[K_LSHIFT]:
                  #pg.image.save(screen,path)
                  pass
                else: #else, save the realscreen
                  pg.image.save(realscreen,path)
                print 'writing to',path
                break #Break loop
                  
          if event.key == K_LEFT: #Move left
            pos = (npos+pos-1)%npos
            loadimage(id,pos)
          elif event.key == K_RIGHT: #Move right
            pos = (pos+1)%npos
            loadimage(id,pos)
          elif event.key == K_UP: #Move up
            id = (id+1)%maxnimpos
            loadimage(id,pos)
          elif event.key == K_DOWN: #Move down
            id = (maxnimpos+id-1)%maxnimpos
            loadimage(id,pos)

          # Numbers
          elif event.key == K_0: 
            goto += '0'
            loadimage(id,pos)
          elif event.key == K_1: 
            goto += '1'
            loadimage(id,pos)
          elif event.key == K_2: 
            goto += '2'
            loadimage(id,pos)
          elif event.key == K_3: 
            goto += '3'
            loadimage(id,pos)
          elif event.key == K_4: 
            goto += '4'
            loadimage(id,pos)
          elif event.key == K_5: 
            goto += '5'
            loadimage(id,pos)
          elif event.key == K_6: 
            goto += '6'
            loadimage(id,pos)
          elif event.key == K_7: 
            goto += '7'
            loadimage(id,pos)
          elif event.key == K_8: 
            goto += '8'
            loadimage(id,pos)
          elif event.key == K_9: 
            goto += '9'
            loadimage(id,pos)
          elif event.key == K_RETURN:
            try:
              index = int(goto)
              if index < maxnimpos:
                id = index
                loadimage(id,pos)
              else:
                print('index out of bounds')
              goto = ''
            except:
              print('index not a number')
              goto = ''
          elif event.key == K_z: 
            # Toggle zoom in zoom out
            pass

      
      screentext = goto if goto=='' else ('load '+goto)
      plot_text(font, screentext, image, (40,30), white)
      realscreen.blit(image, (0, 0))
  
      pg.display.flip() # Flip the display, updates what the user sees
      clock.tick(fps)   # Tick the clock to control fps
  
  while True:
      game(event)
  
