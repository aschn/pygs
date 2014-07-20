#!/usr/bin/env python

import granadata
import os
import wx
import wx.lib.colourselect as csel
import wx.lib.scrolledpanel as scrolled
from math import *


length = 27
radius = 6

wildcard = "text files (*.txt)|*.txt|"\
           "PNG files (*.png)|*.png|"\
           "All files (*.*)|*.*"

class DrawingPanel(scrolled.ScrolledPanel):
    def __init__(self, parent, id):
        scrolled.ScrolledPanel.__init__(self, parent, id)
        self._parent = parent

        self.Bind(wx.EVT_PAINT, self.OnPaint)
        self.Bind(wx.EVT_SET_FOCUS, self.OnSetFocus)
        self.Bind(wx.EVT_LEFT_DOWN, self.OnClick)

        self.SetBackgroundColour(wx.Colour(255, 255, 255))
        self.spacing = 50

        self.SetupScrolling(scroll_y = True)

    def getColor(self, i_type, energy, particle, time):
        pycolor, do_draw_particle, pshape = self._parent.data.getColor(i_type, energy, particle, time)
        if pycolor:
            color = self.tuple2wxColor(pycolor)
        else:
            print "didn't get color for particle", i_type, energy, particle._id
            color = wx.Colour(0,0,0)
        return color, do_draw_particle, pshape

    def tuple2wxColor(self, pycolor):
        if pycolor[0] > 1 or pycolor[1] > 1 or pycolor[2] > 1:
            print "got bad color", pycolor
            return wx.Colour(0,0,0)
        else:
            return wx.Colour(pycolor[0]*255.0, pycolor[1]*255.0, pycolor[2]*255.0)

    def wx2tupleColor(self, wxcolor):
        return (float(wxcolor.Red())/255.0,
                float(wxcolor.Green())/255.0,
                float(wxcolor.Blue())/255.0)

    def OnPaint(self, evt):
        # set up
        dc = wx.PaintDC(self)
        try:
            gc = wx.GraphicsContext.Create(dc)
        except NotImplementedError:
            dc.DrawText("This build of wxPython does not support wx.GraphicsContext", 25, 25)
            return

        font = wx.SystemSettings.GetFont(wx.SYS_DEFAULT_GUI_FONT)
        gc.SetFont(font)
        zoom = self._parent.movieControlPanel.zoom
        gc.Scale(zoom, zoom)
        spacing = self.spacing
        gc.Translate(spacing, spacing)

        time = self._parent.data.params.time
        
        # decide whether to draw overlay
        do_overlay = (2 in self._parent.data.params.layers_to_draw)
        
        # draw boundary rectangles
        layer_width = self._parent.data.params.getWidth(time)
        layer_height = self._parent.data.params.getHeight(time)
        mid_layer_width = 0.5 * layer_width
        mid_layer_height = 0.5 * layer_height
        grana_rad = self._parent.data.params.g_rad
        stroma_width = self._parent.data.params.s_width 
        gc.SetPen(wx.Pen("black"))
        if self._parent.data.params.doEMColors:
            gc.SetBrush(wx.Brush("gray"))
        else:
            gc.SetBrush(wx.Brush("white"))
        for ilayer in range(len(self._parent.data.params.layers_to_draw)):
            gc.DrawRectangle(ilayer*(layer_width+spacing),0, layer_width, layer_height)
        if stroma_width > 0 and grana_rad > 0:
            for ilayer in range(len(self._parent.data.params.layers_to_draw)):
                # draw stroma rectangles
                gc.DrawRectangle(ilayer*(layer_width+spacing),mid_layer_height-0.5*stroma_width, layer_width, stroma_width)
                gc.DrawRectangle(ilayer*(layer_width+spacing)+mid_layer_width-0.5*stroma_width, 0,stroma_width, layer_height)

                # draw grana circles
                gc.DrawEllipse(ilayer*(layer_width+spacing)+mid_layer_width-grana_rad, mid_layer_height-grana_rad, 2.0*grana_rad, 2.0*grana_rad)

        # draw frame
        if time in self._parent.data.trajectory:
            for i_type in range(len(self._parent.data.trajectory[time])):
                type_list = self._parent.data.trajectory[time][i_type]
                gc.PushState()
                for i_layer in range(len(type_list)):
                    layer = type_list[i_layer]
                    is_layer_to_draw = (i_layer in self._parent.data.params.layers_to_draw)
                    if is_layer_to_draw:
                        i_layer_for_drawing = self._parent.data.params.layers_to_draw.index(i_layer)
                    else:
                        i_layer_for_drawing = 0
                    gc.Translate(i_layer_for_drawing * (layer_width+spacing), 0)
                    for i_particle in range(len(layer)):
                        # get position from trajectory
                        particle_pos = layer[i_particle]
                        if particle_pos:
                            x = particle_pos.x
                            y = particle_pos.y
                            particle = self._parent.data.particles[i_type][i_layer][i_particle]
                            energy = particle_pos.energy

                            # make boolean for whether we're going to draw this particle
                            color, do_draw_particle, (theta, pradius, plength) = self.getColor(i_type, energy, particle, time)

                            # clean up particle shape and status
                            theta *= -1.0
                            plength += 2 * pradius
                            tagged = particle._tagged
                            phos = particle._phos

                            # draw particle in own layer
                            if do_draw_particle and is_layer_to_draw:
                                gc.PushState()
                                gc.Translate(x, layer_height - y)
                                gc.Rotate(theta)
                                gc.SetPen(wx.Pen("white", style=wx.TRANSPARENT))
                                gc.SetBrush(wx.Brush(color))
                                gc.DrawRoundedRectangle(-plength/2,-pradius,plength,pradius*2, pradius)
                                gc.PopState()

                                # draw pbc outlines for rods near edges
                                if (plength > 0 and self._parent.data.params.doPBC):
                                    near_edge = False
                                    gc.PushState()
                                    if (x < plength/2): # left edge
                                        near_edge = True
                                        gc.Translate(x + layer_width, layer_height - y)
                                    elif (layer_width - x < plength/2): # right edge
                                        near_edge = True
                                        gc.Translate(x - layer_width, layer_height - y)
                                    elif (y < plength/2): # bottom edge
                                        near_edge = True
                                        gc.Translate(x, - y)
                                    elif (layer_height - y < plength/2): # top edge
                                        near_edge = True
                                        gc.Translate(x, 2*layer_height - y)
                                    if near_edge:
                                        gc.Rotate(theta)
                                        gc.SetPen(wx.Pen(color))
                                        gc.SetBrush(wx.Brush("white", style=wx.TRANSPARENT))
                                        gc.DrawRoundedRectangle(-plength/2,-pradius,plength,pradius*2, pradius)
                                    gc.PopState()

                                # draw PSII cores
                                if i_type == 1 and self._parent.data.params.doEMColors:
                                    clength = 15.3
                                    cradius = 3.8
                                    ctheta = theta - (2*i_layer - 1) * 0.389
                                    gc.PushState()
                                    gc.Translate(x, layer_height - y)
                                    gc.Rotate(ctheta)
                                    gc.SetPen(wx.Pen("black"))
                                    gc.SetBrush(wx.Brush("white", style=wx.TRANSPARENT))
                                    gc.DrawRoundedRectangle(-clength/2,-cradius,clength,cradius*2, cradius)
                                    gc.PopState()
                            
                                            
                            # draw particle in overlay
                            if do_draw_particle and do_overlay:
                                gc.PushState()
                                gc.Translate((len(self._parent.data.params.layers_to_draw) - 1 - i_layer_for_drawing) * (layer_width+spacing), 0)
                                gc.Translate(x, layer_height - y)
                                gc.Rotate(theta)
                                gc.SetPen(wx.Pen(color))
                                gc.SetBrush(wx.Brush("white", style=wx.TRANSPARENT))
                                gc.DrawRoundedRectangle(-plength/2,-pradius,plength,pradius*2, pradius)
                                gc.PopState()
                gc.PopState()
                            
        else:
            print "can't draw time", time
            
    def OnSetFocus(self, event):
        self.Refresh()

    def OnClick(self, event):
        # get frame info
        time = self._parent.data.params.time
        zoom = self._parent.movieControlPanel.zoom
        layer_width = self._parent.data.params.getWidth(time)
        layer_height = self._parent.data.params.getHeight(time)
        spacing = self.spacing

        # get mouse info
        (mx, my) = event.GetPositionTuple()
        mx /= zoom
        my /= zoom

        # get layer
        layer = -1
        if my > spacing and my < spacing + layer_height: # valid y
            my = layer_height + spacing - my # translate y to simulation coords
            if mx > spacing and mx < spacing + layer_width:
                layer = 0
                mx = mx - spacing # translate x to simulation coords
                print "looking in layer 0"
            elif mx > spacing * 2 + layer_width and mx < (spacing + layer_width) * 2:
                layer = 1
                mx = mx - layer_width - 2 * spacing # translate x to simulation coords
                print "looking in layer 1"
            else:
                print mx, my, "not in valid x region"
                return False
        else:
            print mx, my, "not in valid y region"
            return False

        if layer > -1:
            # check free LHCIIs
            lhccoords = self._parent.data.trajectory[time][0][layer]
            for lhc in self._parent.data.particles[0][layer]:
                if lhc:
                    if lhccoords[lhc._id]: # skip if missing that id
                        px = lhccoords[lhc._id].x
                        py = lhccoords[lhc._id].y
                        if (px-mx)**2 + (py-my)**2 < lhc._radius**2:
                            lhc._tagged = not lhc._tagged 
                            self.Refresh()
                            print "toggled free lhc", lhc._id, "at", px, py
                            return True

            # check PSIIs
            pscoords = self._parent.data.trajectory[time][1][layer]
            for ps in self._parent.data.particles[1][layer]:
                if ps:
                    if pscoords[ps._id]: # skip if missing that id
                        px = pscoords[ps._id].x
                        py = pscoords[ps._id].y
                        if (px-mx)**2 + (py-my)**2 < (ps._length / 2)**2:
                            ps._tagged = not ps._tagged 
                            self.Refresh()
                            print "toggled ps", ps._id, "at", px, py
                            return True

            # check bound LHCIIs
            lhccoords = self._parent.data.trajectory[time][2][layer]
            for lhc in self._parent.data.particles[2][layer]:
                if lhc:
                    if lhccoords[lhc._id]: # skip if missing that id
                        px = lhccoords[lhc._id].x
                        py = lhccoords[lhc._id].y
                        if (px-mx)**2 + (py-my)**2 < lhc._radius**2:
                            lhc._tagged = not lhc._tagged 
                            self.Refresh()
                            print "toggled bound lhc", lhc._id
                            return True

            # check discs
            disccoords = self._parent.data.trajectory[time][5][layer]
            for disc in self._parent.data.particles[5][layer]:
                if disc: # skip if missing that id
                    if disccoords[disc._id]:
                        px = disccoords[disc._id].x
                        py = disccoords[disc._id].y
                        if (px-mx)**2 + (py-my)**2 < disc._radius**2:
                            disc._tagged = not disc._tagged 
                            self.Refresh()
                            print "toggled bound disc", disc._id
                            return True

class ColorControlPanel(wx.Panel):
    def __init__(self, parent, id, painter):
        self._parent = parent
        self._painter = painter
        wx.Panel.__init__(self, parent, id)
        self.SetAutoLayout(True)
        mainSizer = wx.FlexGridSizer(4,2)
        self.SetSizer(mainSizer)

        colorbox = wx.StaticBox(self, -1, 'colors')
        colorBoxSizer = wx.StaticBoxSizer(colorbox)
        colorSizer = wx.FlexGridSizer(1, 2) # sizer to contain all the example buttons
       
        # build several examples of buttons with different colours and sizes
        buttonData = [
            (self.OnSelectColorLHCDefault,
             self._parent.drawingPanel.tuple2wxColor(self._parent.data.params.lhcDefaultColor),
             wx.DefaultSize, "LHCII"),
            (self.OnSelectColorLHCTag,
             self._parent.drawingPanel.tuple2wxColor(self._parent.data.params.lhcTagColor),
             wx.DefaultSize, "LHCII tag"),
            (self.OnSelectColorPSDefault,
             self._parent.drawingPanel.tuple2wxColor(self._parent.data.params.psDefaultColor),
             wx.DefaultSize, "PSII"),
            (self.OnSelectColorPSTag,
             self._parent.drawingPanel.tuple2wxColor(self._parent.data.params.psTagColor),
             wx.DefaultSize, "PSII tag"),
            ]

        self.buttonRefs = [] # for saving references to buttons

        # build each button and save a reference to it
        for event, color, size, label in buttonData:
            b = csel.ColourSelect(self, -1, label, color, size = size)
            b.Bind(csel.EVT_COLOURSELECT, event)
            self.buttonRefs.append(("", b))  # store reference to button
            colorSizer.AddMany([
                (b, 0, wx.ALL, 3),
                ])
        colorBoxSizer.Add(colorSizer)
        mainSizer.Add(colorBoxSizer, 0, wx.ALL, 3)

        # build color scheme chooser radiobox
        schemebox = wx.StaticBox(self, -1, 'color schemes')
        schemeBoxSizer = wx.StaticBoxSizer(schemebox)
        schemeStaticSizer = wx.FlexGridSizer(4,2)
        self.tagToggle = wx.ToggleButton(self, -1, 'use tags')
        self.tagToggle.SetValue(self._parent.data.params.doTagColors)
        self.Bind(wx.EVT_TOGGLEBUTTON, self.OnTagToggle, self.tagToggle)
        schemeStaticSizer.Add(self.tagToggle)

        self.trimerEnergyToggle = wx.ToggleButton(self, -1, 'tri. energy')
        self.trimerEnergyToggle.SetValue(self._parent.data.params.doTrimerEnergyColors)
        self.Bind(wx.EVT_TOGGLEBUTTON, self.OnTrimerEnergyToggle, self.trimerEnergyToggle)
        schemeStaticSizer.Add(self.trimerEnergyToggle)

   #     self.psStackToggle = wx.ToggleButton(self, -1, 'RC energy')
   #     self.psStackToggle.SetValue(self._parent.data.params.doPSStackColors)
   #     self.Bind(wx.EVT_TOGGLEBUTTON, self.OnPSStackToggle, self.psStackToggle)
   #     schemeStaticSizer.Add(self.psStackToggle)

        self.mTrimerToggle = wx.ToggleButton(self, -1, 'M-LHCII')
        self.mTrimerToggle.SetValue(self._parent.data.params.doMTrimerColors)
        self.Bind(wx.EVT_TOGGLEBUTTON, self.OnMTrimerToggle, self.mTrimerToggle)
        schemeStaticSizer.Add(self.mTrimerToggle)

    #    self.monomerEnergyToggle = wx.ToggleButton(self, -1, 'mon. energy')
    #    self.monomerEnergyToggle.SetValue(self._parent.data.params.doMonomerEnergyColors)
    #    self.Bind(wx.EVT_TOGGLEBUTTON, self.OnMonomerEnergyToggle, self.monomerEnergyToggle)
    #    schemeStaticSizer.Add(self.monomerEnergyToggle)

        self.xtalToggle = wx.ToggleButton(self, -1, 'xtal')
        self.xtalToggle.SetValue(self._parent.data.params.doXtalColors)
        self.Bind(wx.EVT_TOGGLEBUTTON, self.OnXtalToggle, self.xtalToggle)
        schemeStaticSizer.Add(self.xtalToggle)

        self.shapeToggle = wx.ToggleButton(self, -1, 'shape')
        self.shapeToggle.SetValue(self._parent.data.params.doDelaunayColors)
        self.Bind(wx.EVT_TOGGLEBUTTON, self.OnShapeToggle, self.shapeToggle)
        schemeStaticSizer.Add(self.shapeToggle)

     #   self.dispToggle = wx.ToggleButton(self, -1, 'disp')
     #   self.dispToggle.SetValue(self._parent.data.params.doDispColors)
     #   self.Bind(wx.EVT_TOGGLEBUTTON, self.OnDispToggle, self.dispToggle)
     #   schemeStaticSizer.Add(self.dispToggle)

        self.megaEnergyToggle = wx.ToggleButton(self, -1, 'mega. energy')
        self.megaEnergyToggle.SetValue(self._parent.data.params.doMsiteEnergyColors)
        self.Bind(wx.EVT_TOGGLEBUTTON, self.OnMegaEnergyToggle, self.megaEnergyToggle)
        schemeStaticSizer.Add(self.megaEnergyToggle)

   #     self.EMToggle = wx.ToggleButton(self, -1, 'EM')
   #     self.EMToggle.SetValue(self._parent.data.params.doEMColors)
   #     self.Bind(wx.EVT_TOGGLEBUTTON, self.OnEMToggle, self.EMToggle)
   #     schemeStaticSizer.Add(self.EMToggle)
        
        self.HexToggle = wx.ToggleButton(self, -1, 'hexatic')
        self.HexToggle.SetValue(self._parent.data.params.doHexaticColors)
        self.Bind(wx.EVT_TOGGLEBUTTON, self.OnHexToggle, self.HexToggle)
        schemeStaticSizer.Add(self.HexToggle)

        self.clusterToggle = wx.ToggleButton(self, -1, 'cluster')
        self.clusterToggle.SetValue(self._parent.data.params.doClusterColors)
        self.Bind(wx.EVT_TOGGLEBUTTON, self.OnClusterToggle, self.clusterToggle)
        schemeStaticSizer.Add(self.clusterToggle)

        schemeBoxSizer.Add(schemeStaticSizer)
        mainSizer.Add(schemeBoxSizer)
#        self.schemenames = ["energy", "tags"]
#        self.colorSchemeRadio = wx.RadioBox(
#                self, -1, "color scheme", wx.DefaultPosition, wx.DefaultSize,
#                self.schemenames, 1, wx.RA_SPECIFY_ROWS
#                )
#        
#        self.Bind(wx.EVT_RADIOBOX, self.OnRadioBox, self.colorSchemeRadio)
#        mainSizer.Add(self.colorSchemeRadio)

        # build "tag rectangle" buttons
        rectbox = wx.StaticBox(self, -1, 'tag within a rectangle')
        rectBoxSizer = wx.StaticBoxSizer(rectbox)
        tagRectSizer = wx.FlexGridSizer(3, 2)
        self.xmin = wx.TextCtrl(self, -1, "xmin", (-1, -1), (50, -1) )
        tagRectSizer.Add(self.xmin)
        self.xmax = wx.TextCtrl(self, -1, "xmax", (-1, -1), (50, -1) )
        tagRectSizer.Add(self.xmax)
        self.ymin = wx.TextCtrl(self, -1, "ymin", (-1, -1), (50, -1) )
        tagRectSizer.Add(self.ymin)
        self.ymax = wx.TextCtrl(self, -1, "ymax", (-1, -1), (50, -1) )
        tagRectSizer.Add(self.ymax)
        tagRect0 = wx.Button(self, -1, "bottom")
        self.Bind(wx.EVT_BUTTON, self.OnClickTagRect0, tagRect0)
        tagRectSizer.Add(tagRect0)
        tagRect1 = wx.Button(self, -1, "top")
        self.Bind(wx.EVT_BUTTON, self.OnClickTagRect1, tagRect1)
        tagRectSizer.Add(tagRect1)
        rectBoxSizer.Add(tagRectSizer)
        mainSizer.Add(rectBoxSizer, 0, wx.ALL, 3)

        # build "tag circle" buttons
        circlebox = wx.StaticBox(self, -1, 'tag within a circle')
        circleBoxSizer = wx.StaticBoxSizer(circlebox, wx.VERTICAL)
        tagCircleSizer = wx.FlexGridSizer(3,2)
        self.circlex = wx.TextCtrl(self, -1, "x", (-1, -1), (50, -1) )
        tagCircleSizer.Add(self.circlex)
        self.circley = wx.TextCtrl(self, -1, "y", (-1, -1), (50, -1) )
        tagCircleSizer.Add(self.circley)
        self.circlerad = wx.TextCtrl(self, -1, "radius", (-1, -1), (50, -1) )
        tagCircleSizer.Add(self.circlerad)
        tagCircle0 = wx.Button(self, -1, "bottom")
        self.Bind(wx.EVT_BUTTON, self.OnClickTagCircle0, tagCircle0)
        tagCircleSizer.Add(tagCircle0)
        tagCircle1 = wx.Button(self, -1, "top")
        self.Bind(wx.EVT_BUTTON, self.OnClickTagCircle1, tagCircle1)
        tagCircleSizer.Add(tagCircle1)
        circleBoxSizer.Add(tagCircleSizer)
        mainSizer.Add(circleBoxSizer, 0, wx.ALL, 3)

        # build "untag all" buttons
        untagbox = wx.StaticBox(self, -1, 'untag all')
        untagSizer = wx.StaticBoxSizer(untagbox)
        untag0 = wx.Button(self, -1, "bottom")
        self.Bind(wx.EVT_BUTTON, self.OnClickUntag0, untag0)
        untagSizer.Add(untag0)
        untag1 = wx.Button(self, -1, "top")
        self.Bind(wx.EVT_BUTTON, self.OnClickUntag1, untag1)
        untagSizer.Add(untag1)
        mainSizer.Add(untagSizer, 0, wx.ALL, 3)

        self.Layout()

    def OnSelectColorLHCDefault(self, event):
         self._parent.data.params.lhcDefaultColor = self._parent.drawingPanel.wx2tupleColor(event.GetValue())
         self._painter.Refresh(False)
    def OnSelectColorLHCTag(self, event):
         self._parent.data.params.lhcTagColor = self._parent.drawingPanel.wx2tupleColor(event.GetValue())
         self._painter.Refresh(False)
    def OnSelectColorPSDefault(self, event):
         self._parent.data.params.psDefaultColor = self._parent.drawingPanel.wx2tupleColor(event.GetValue())
         self._painter.Refresh(False)
    def OnSelectColorPSTag(self, event):
         self._parent.data.params.psTagColor = self._parent.drawingPanel.wx2tupleColor(event.GetValue())
         self._painter.Refresh(False)

    def OnRadioBox(self, event):
        i_radio = self.colorSchemeRadio.GetSelection()
        if i_radio==0:
            self._parent.data.params.doEnergyColors = True
            self._parent.data.params.doTagColors = False
        elif i_radio==1:
            self._parent.data.params.doEnergyColors = False
            self._parent.data.params.doTagColors = True
        else:
            print "got invalid selection from color scheme selector"
        self._painter.Refresh()

    def OnTagToggle(self, event):
        self._parent.data.params.doTagColors = not self._parent.data.params.doTagColors
        self._painter.Refresh()

    def OnTrimerEnergyToggle(self, event):
        self._parent.data.params.doTrimerEnergyColors = not self._parent.data.params.doTrimerEnergyColors
        self._painter.Refresh()

    def OnMonomerEnergyToggle(self, event):
        self._parent.data.params.doMonomerEnergyColors = not self._parent.data.params.doMonomerEnergyColors
        self._painter.Refresh()

    def OnMegaEnergyToggle(self, event):
        self._parent.data.params.doMsiteEnergyColors = not self._parent.data.params.doMsiteEnergyColors
        self._painter.Refresh()

    def OnPSStackToggle(self, event):
        self._parent.data.params.doPSStackColors = not self._parent.data.params.doPSStackColors
        self._painter.Refresh()

    def OnXtalToggle(self, event):
        self._parent.data.params.doXtalColors = not self._parent.data.params.doXtalColors
        self._painter.Refresh()

    def OnDispToggle(self, event):
        self._parent.data.params.doDispColors = not self._parent.data.params.doDispColors
        self._painter.Refresh()
        
    def OnEMToggle(self, event):
        self._parent.data.params.doEMColors = not self._parent.data.params.doEMColors
        self._painter.Refresh()

    def OnClusterToggle(self, event):
        self._parent.data.params.doClusterColors = not self._parent.data.params.doClusterColors
        self._painter.Refresh()

    def OnShapeToggle(self, event):
        self._parent.data.params.doDelaunayColors = not self._parent.data.params.doDelaunayColors
        self._painter.Refresh()

    def OnHexToggle(self, event):
        self._parent.data.params.doHexaticColors = not self._parent.data.params.doHexaticColors
        self._painter.Refresh()

    def OnMTrimerToggle(self, event):
        self._parent.data.params.doMTrimerColors = not self._parent.data.params.doMTrimerColors
        self._painter.Refresh()

    def OnClickTagRect0(self, event):
        self.tagRect(0)
    def OnClickTagRect1(self, event):
        self.tagRect(1)

    def OnClickTagCircle0(self, event):
        self.tagCircle(0)
    def OnClickTagCircle1(self, event):
        self.tagCircle(1)

    def OnClickUntag0(self, event):
        self.untag(0)
    def OnClickUntag1(self, event):
        self.untag(1)

    def tagRect(self, layer):
        time = self._parent.data.params.time
        xmin = float(self.xmin.GetValue())
        xmax = float(self.xmax.GetValue())
        ymin = float(self.ymin.GetValue())
        ymax = float(self.ymax.GetValue())
        for type in (0, 1, 2):
            for particle in self._parent.data.particles[type][layer]:
                if particle:
                    if self._parent.data.trajectory[time][type][layer][particle._id]:
                        point = self._parent.data.trajectory[time][type][layer][particle._id]
                        if point.x > xmin and point.x < xmax and point.y > ymin and point.y < ymax:
                            particle._tagged = True
        self._painter.Refresh()

    def tagCircle(self, layer):
        time = self._parent.data.params.time
        x = float(self.circlex.GetValue())
        y = float(self.circley.GetValue())
        rad = float(self.circlerad.GetValue())
        for type in (0, 1, 2):
            for particle in self._parent.data.particles[type][layer]:
                if particle:
                    if self._parent.data.trajectory[time][type][layer][particle._id]:
                        point = self._parent.data.trajectory[time][type][layer][particle._id]
                        if (point.x - x)**2 + (point.y  - y)**2 < rad**2:
                            particle._tagged = True
        self._painter.Refresh()

    def untag(self, layer):
        for type in (0, 1, 2):
            for particle in self._parent.data.particles[type][layer]:
                if particle:
                    particle._tagged = False
        self._painter.Refresh()



class MovieControlPanel(wx.Panel):
    def __init__(self, parent, id):
        wx.Panel.__init__(self, parent, id)
        self._parent = parent
        vbox = wx.FlexGridSizer(rows=4, cols=1, vgap=10)

        # build timer
        self.timer = wx.Timer(self, -1)
        self.Bind(wx.EVT_TIMER, self.OnTimer, self.timer)
        self.fps = 100

        # set initial time
        #self.time = 0
        #self.timechoices = []
        #self.stride = 1

        # build play/pause button
        self.playButton = wx.Button(self, -1, "play")
        self.Bind(wx.EVT_BUTTON, self.OnClickPlay, self.playButton)
        vbox.Add(self.playButton, 1, wx.EXPAND) 

        # build draw frame button
        psSizer = wx.FlexGridSizer(1,2)
        self.psDrawButton = wx.Button(self, -1, "draw PS frame")
        self.Bind(wx.EVT_BUTTON, self.OnClickDraw, self.psDrawButton)
        psSizer.Add(self.psDrawButton, 1, wx.EXPAND) 
        self.psMovieToggle = wx.ToggleButton(self, -1, 'record PS movie')
        self.psMovieToggle.SetValue(False)
        self.Bind(wx.EVT_TOGGLEBUTTON, self.OnMovieToggle, self.psMovieToggle)
        psSizer.Add(self.psMovieToggle)        
        vbox.Add(psSizer, 1, wx.EXPAND)

        # build layer toggles
        layerbox = wx.StaticBox(self, -1, 'layers')
        layerBoxSizer = wx.StaticBoxSizer(layerbox)
        self.botToggle = wx.ToggleButton(self, -1, 'bottom')
        self.botToggle.SetValue(True)
        self.Bind(wx.EVT_TOGGLEBUTTON, self.OnLayerToggle, self.botToggle)
        layerBoxSizer.Add(self.botToggle)        
        self.topToggle = wx.ToggleButton(self, -1, 'top')
        self.topToggle.SetValue(True)
        self.Bind(wx.EVT_TOGGLEBUTTON, self.OnLayerToggle, self.topToggle)
        layerBoxSizer.Add(self.topToggle)
        self.olToggle = wx.ToggleButton(self, -1, 'overlay')
        self.olToggle.SetValue(True)
        self.Bind(wx.EVT_TOGGLEBUTTON, self.OnLayerToggle, self.olToggle)
        layerBoxSizer.Add(self.olToggle)                
        vbox.Add(layerBoxSizer, 1, wx.EXPAND)
       
        # build time info
        timestaticbox = wx.StaticBox(self, -1, "movie controls")
        timeStaticSizer = wx.StaticBoxSizer(timestaticbox)
        timebox = wx.FlexGridSizer(rows=5, cols=2, hgap=10, vgap=5)
        timebox.AddGrowableCol(1)
        # min time
        timebox.Add(wx.StaticText(self, -1, "min time (us)"), wx.CENTER)
        self.timeMin = wx.StaticText(self, -1, "100000")
        timebox.Add(self.timeMin, wx.CENTER)
        # max time
        timebox.Add(wx.StaticText(self, -1, "max time (us)"), wx.CENTER)
        self.timeMax = wx.StaticText(self, -1, "-100000")
        timebox.Add(self.timeMax, wx.CENTER)
        # stride
        timebox.Add(wx.StaticText(self, -1, "stride (us)"), wx.CENTER)
        self.timeStride = wx.TextCtrl(self, -1, str(self._parent.data.params.stride))
        self.Bind(wx.EVT_TEXT, self.OnTimeStrideUpdate, self.timeStride)
        timebox.Add(self.timeStride, wx.CENTER)
        # fps
#        timebox.Add(wx.StaticText(self, -1, "nominal fps"), wx.CENTER)
#        self.timeFps = wx.TextCtrl(self, -1, str(self.fps))
#        timebox.Add(self.timeFps, wx.CENTER)
        # current time
        timebox.Add(wx.StaticText(self, -1, "current time"), wx.CENTER)
        strtimes = [str(time) for time in self._parent.data.params.times]
        self.timeCombo = wx.ComboBox(self, -1, '0',(-1,-1), (-1,-1), strtimes, wx.CB_SORT)
        self.Bind(wx.EVT_TEXT, self.OnTimeComboUpdate, self.timeCombo)
        timebox.Add(self.timeCombo, 1, wx.EXPAND)

        # build zoom combobox in timebox
        self.zoom = 0.75
        timebox.Add(wx.StaticText(self, -1, "zoom"), wx.CENTER)
        zoomchoices = ['0.1', '0.25', '0.5', '0.75', '1.0', '1.5', '2.0']
        self.zoomCombo = wx.ComboBox(self, -1, '0.75',(-1,-1), (100,-1), zoomchoices)
        self.Bind(wx.EVT_TEXT, self.OnZoom, self.zoomCombo)
        timebox.Add(self.zoomCombo, 1, wx.EXPAND)

        # finish adding to sizers
        timeStaticSizer.Add(timebox)
        vbox.Add(timeStaticSizer, 1, wx.EXPAND)

        # build redraw button
#        redrawButton = wx.Button(self, -1, "redraw")
#        self.Bind(wx.EVT_BUTTON, self.OnClickRedraw, redrawButton)
#        vbox.Add(redrawButton, 1, wx.EXPAND)

        # add sizer to panel
        self.SetSizer(vbox)

    def OnClickPlay(self, event):
        btnLabel = self.playButton.GetLabel()
        if btnLabel == "play":
            print "starting movie..."
            self.timer.Start(1000 / int(self.fps)) # 1000 / fps = milliseconds per frame
            self.playButton.SetLabel("pause")
        else:
            print "movie stopped!"
            self.timer.Stop()
            self.playButton.SetLabel("play")

    def OnClickDraw(self, event):
        self._parent.data.DrawPS(int(self.timeCombo.GetValue()))

    def OnMovieToggle(self, event):
        print "will record frames as postscript files"

    def OnLayerToggle(self, event):
        layers = []
        if self.botToggle.GetValue():
            layers.append(0)
        if self.topToggle.GetValue():
            layers.append(1)
        if self.olToggle.GetValue():
            layers.append(2)
        self._parent.data.params.layers_to_draw = layers
        self._parent.drawingPanel.Refresh()

    def OnTimeComboUpdate(self, event):
        time = int(self.timeCombo.GetValue())
        self.SetTime(time)

    def OnTimeStrideUpdate(self, event):
        self._parent.data.params.setStride(int(self.timeStride.GetValue()))
        self._parent.drawingPanel.Refresh()

    def SetTime(self, time):
        self._parent.data.params.setTime(int(time))
        self.timeCombo.SetValue(str(time))
        self._parent.fileControlPanel.UpdateXY()
        if self._parent.data.params.doMTrimerColors:
            self._parent.data.TagMtrimers(time)
        self._parent.drawingPanel.Refresh()
        if (self.psMovieToggle.GetValue()):
            if self.time % self.stride == 0:
                self._parent.data.DrawPS(time)

    def OnTimer(self, event):
#        self.timeFps.SetLabel(str(1000 / self.timer.GetInterval()))  # 1000 / fps = milliseconds per framexs
        if self._parent.data.params.time in self._parent.data.params.times:
            i_time = self._parent.data.params.times.index(self._parent.data.params.time)
            if i_time < len(self._parent.data.params.times) - 1 and self._parent.data.params.times[i_time + 1] - self._parent.data.params.time == self._parent.data.params.stride:
                self.SetTime(self._parent.data.params.times[i_time + 1])
            else:
                print "end of movie"
                self.timer.Stop()
                self.playButton.SetLabel("play")

        else:
            print "not at valid time, can't play movie"
            self.timer.Stop()
            self.playButton.SetLabel("play")

    def UpdateTimechoices(self):
        for timechoice in self._parent.data.params.times:
            if timechoice != int(self.timeCombo.GetValue()):
                self.timeCombo.Append(str(timechoice))
            if timechoice > int(self.timeMax.GetLabel()):
                self.timeMax.SetLabel(str(timechoice))
            elif timechoice < int(self.timeMin.GetLabel()):
                self.timeMin.SetLabel(str(timechoice))
        self._parent.data.params.setStrideAuto()
        self.timeStride.SetValue(str(self._parent.data.params.stride))

    def OnZoom(self, event):
        self.zoom = float(self.zoomCombo.GetValue())
        self._parent.drawingPanel.Refresh()
        
class FileControlPanel(wx.Panel):
    def __init__(self, parent, id):
        wx.Panel.__init__(self, parent, id)
        self._parent = parent

    #    self.width = 266.0
    #    self.height = 266.0
    #    self.g_rad = 0.0
    #    self.s_width = 0.0
        self.epsilon = 1.0
        self.mega = 0.0

        staticbox = wx.StaticBox(self, -1, 'initialize data')
        boxsizer = wx.StaticBoxSizer(staticbox, wx.VERTICAL)
        sizer = wx.FlexGridSizer(cols=1, hgap=10, vgap=5)

        # build file chooser
        filebutton = wx.Button(self, -1, "open files")
        self.Bind(wx.EVT_BUTTON, self.OnButton, filebutton)
        sizer.Add(filebutton, 1, wx.ALIGN_CENTER)

        # build directory label
        self.dir_text = wx.StaticText(self, -1, '')
        sizer.Add(self.dir_text, 1, wx.ALIGN_LEFT)

        # build x-y and energy text controls
        xysizer = wx.FlexGridSizer(rows=3, cols=4, hgap=10, vgap=5)
        xysizer.Add(wx.StaticText(self, -1, 'width ='), 1, wx.EXPAND)
        self.xcntrl = wx.TextCtrl(self, -1, str(self._parent.data.params.width), (-1, -1), (50, -1) )
        self.Bind(wx.EVT_TEXT, self.OnX, self.xcntrl)
        xysizer.Add(self.xcntrl)

        xysizer.Add(wx.StaticText(self, -1, 'height ='), 1, wx.EXPAND)
        self.ycntrl = wx.TextCtrl(self, -1, str(self._parent.data.params.height), (-1, -1), (50, -1) )
        self.Bind(wx.EVT_TEXT, self.OnY, self.ycntrl)
        xysizer.Add(self.ycntrl)

        xysizer.Add(wx.StaticText(self, -1, 'grana rad ='), 1, wx.EXPAND)
        self.grcntrl = wx.TextCtrl(self, -1, str(self._parent.data.params.g_rad), (-1, -1), (50, -1) )
        self.Bind(wx.EVT_TEXT, self.OnGR, self.grcntrl)
        xysizer.Add(self.grcntrl)

        xysizer.Add(wx.StaticText(self, -1, 'stroma wdth ='), 1, wx.EXPAND)
        self.swcntrl = wx.TextCtrl(self, -1, str(self._parent.data.params.s_width), (-1, -1), (50, -1) )
        self.Bind(wx.EVT_TEXT, self.OnSW, self.swcntrl)
        xysizer.Add(self.swcntrl)


        xysizer.Add(wx.StaticText(self, -1, 'e_stroma ='), 1, wx.EXPAND)
        self.epsiloncntrl = wx.TextCtrl(self, -1, str(self.epsilon), (-1, -1), (50, -1) )
        self.Bind(wx.EVT_TEXT, self.OnEpsilon, self.epsiloncntrl)
        xysizer.Add(self.epsiloncntrl)

        xysizer.Add(wx.StaticText(self, -1, 'e_mega ='), 1, wx.EXPAND)
        self.megacntrl = wx.TextCtrl(self, -1, str(self.mega), (-1, -1), (50, -1) )
        self.Bind(wx.EVT_TEXT, self.OnMega, self.megacntrl)
        xysizer.Add(self.megacntrl)

        sizer.Add(xysizer, 1, wx.EXPAND)

        # build PBC toggle
        self.pbcToggle = wx.ToggleButton(self, -1, 'use PBC')
        self.pbcToggle.SetValue(self._parent.data.params.doPBC)
        self.Bind(wx.EVT_TOGGLEBUTTON, self.OnPbcToggle, self.pbcToggle)
        sizer.Add(self.pbcToggle, 1, wx.ALIGN_CENTER)

        boxsizer.Add(sizer)
        self.SetSizer(boxsizer)

    def OnButton(self, event):
        # Create the dialog. In this case the current directory is forced as the starting
        # directory for the dialog, and no default file name is forced. This can easilly
        # be changed in your program. This is an 'open' dialog, and allows multitple
        # file selections as well.
        #
        # Finally, if the directory is changed in the process of getting files, this
        # dialog is set up to change the current working directory to the path chosen.
        dlg = wx.FileDialog(
            self, message="Choose a file",
            defaultDir=os.getcwd(), 
            defaultFile="config*.txt",
            wildcard=wildcard,
            style=wx.OPEN | wx.MULTIPLE | wx.CHANGE_DIR
            )

        # Show the dialog and retrieve the user response. If it is the OK response, 
        # process the data.
        if dlg.ShowModal() == wx.ID_OK:
            # This returns a Python list of files that were selected.
            paths = dlg.GetPaths()

            for path in paths:
                self.readFile(path)

            # update wrt sample path
            full_path = os.path.dirname(paths.pop())
            self._parent.data.params.setOutputPath(full_path)
            sampledirs = full_path.split('/')
            last_one = sampledirs.pop()
            last_two = sampledirs.pop() + '/' + last_one
            self.dir_text.SetLabel(last_two)
            self._parent.movieControlPanel.UpdateTimechoices()
            

        # Destroy the dialog. Don't do this until you are done with it!
        # BAD things can happen otherwise!
        dlg.Destroy()

    def readFile(self, filename):
        self._parent.data.readFile(filename)

    def OnX(self, event):
        self._parent.data.params.setWidth(float(self.xcntrl.GetValue()))
        self._parent.drawingPanel.Refresh()
  #      self.UpdateXY()

    def OnY(self, event):
        self._parent.data.params.setHeight(float(self.ycntrl.GetValue()))
        self._parent.drawingPanel.Refresh()
  #      self.UpdateXY()

    def OnGR(self, event):
        self._parent.data.params.setGRad(float(self.grcntrl.GetValue()))
        self._parent.drawingPanel.Refresh()

    def OnSW(self, event):
        self._parent.data.params.setSWidth(float(self.swcntrl.GetValue()))
        self._parent.drawingPanel.Refresh()
        
    def OnEpsilon(self, event):
        self._parent.data.params.lhc_stacking_epsilon = float(self.epsiloncntrl.GetValue())
        self._parent.drawingPanel.Refresh()

    def OnMega(self, event):
        self._parent.data.params.megacomplex_epsilon = float(self.megacntrl.GetValue())
        self._parent.drawingPanel.Refresh()

    def UpdateXY(self):
        self.xcntrl.SetValue(str(self._parent.data.params.width))
        self.ycntrl.SetValue(str(self._parent.data.params.height))
#       try:
 #           self.xcntrl.SetValue('%d' % self.width)
 #           self.ycntrl.SetValue('%d' % self.height)
 #           self._parent.drawingPanel.Refresh()
 #       except RuntimeError:
 #           print "runtime error: width", '%d' % self.width, "and height", '%d' % self.height

    def OnPbcToggle(self, event):
        self._parent.data.params.doPBC = not self._parent.data.params.doPBC
        self._parent.drawingPanel.Refresh()

class PlotControlPanel(wx.Panel):
    def __init__(self, parent, id):
        wx.Panel.__init__(self, parent, id)
        self._parent = parent

        # build sizer
        plotbox = wx.StaticBox(self, -1, 'plot correlation fxns')
        plotBoxSizer = wx.StaticBoxSizer(plotbox)
        plotSizer = wx.FlexGridSizer(1, 2)

        # build buttons
        self.msd_lhc = wx.Button(self, -1, "LHCII MSD", (-1, -1), (100, -1) )
        self.Bind(wx.EVT_BUTTON, self.OnMsdLhc, self.msd_lhc)
        plotSizer.Add(self.msd_lhc)
        self.msd_ps = wx.Button(self, -1, "PSII MSD", (-1, -1), (100, -1) )
        self.Bind(wx.EVT_BUTTON, self.OnMsdPs, self.msd_ps)
        plotSizer.Add(self.msd_ps)
        self.nn_ps = wx.Button(self, -1, "PSII NN", (-1, -1), (100, -1) )
        self.Bind(wx.EVT_BUTTON, self.OnNNPs, self.nn_ps)
        plotSizer.Add(self.nn_ps)
        self.regions = wx.Button(self, -1, "regions", (-1, -1), (100, -1) )
        self.Bind(wx.EVT_BUTTON, self.OnRegions, self.regions)
        plotSizer.Add(self.regions)
   #     self.gofr_lhc = wx.Button(self, -1, "LHCII g(r)", (-1, -1), (100, -1) )
   #     self.Bind(wx.EVT_BUTTON, self.OnGofRLhc, self.gofr_lhc)
   #     plotSizer.Add(self.gofr_lhc)
   #     self.gofr_ps = wx.Button(self, -1, "PSII g(r,th)", (-1, -1), (100, -1) )
   #     self.Bind(wx.EVT_BUTTON, self.OnGofTheta, self.gofr_ps)
#        plotSizer.Add(self.gofr_ps)
        self.gofr_ps = wx.Button(self, -1, "PSII g(r)", (-1, -1), (100, -1) )
        self.Bind(wx.EVT_BUTTON, self.OnGofRPs, self.gofr_ps)
        plotSizer.Add(self.gofr_ps)
        self.cluster = wx.Button(self, -1, "clusters", (-1, -1), (100, -1) )
        self.Bind(wx.EVT_BUTTON, self.OnClusters, self.cluster)
        plotSizer.Add(self.cluster)

        plotBoxSizer.Add(plotSizer)
        self.SetSizer(plotBoxSizer)

    def OnMsdLhc(self, event):
        self._parent.data.PlotMsd(int(0))

    def OnMsdPs(self, event):
        self._parent.data.PlotMsd(1)

    def OnNNPs(self, event):
        self._parent.data.PlotNearestNbhrs(1)

    def OnRegions(self, event):
        self._parent.data.PlotRegions()

    def OnGofRLhc(self, event):
        self._parent.data.PlotGofR(0)

    def OnGofRPs(self, event):
        self._parent.data.PlotGofR(1)

    def OnGofTheta(self, event):
        self._parent.data.PlotGofTheta(1)

    def OnClusters(self, event):
        #self.PlotClusters()
        self._parent.data.PlotArrays()


class MainFrame(wx.Frame):
    def __init__(self, parent, id, title):
        wx.Frame.__init__(self, parent, id, title, size=(1000,750))

        # build parent structure, not sure why it needs to be panel not self
        panel = wx.Panel(self, -1)

        # initialize data structure
        panel.data = granadata.GranaData()
        
        # build panels
        panel.drawingPanel = DrawingPanel(panel, -1)
        panel.movieControlPanel = MovieControlPanel(panel, -1)
        panel.colorControlPanel = ColorControlPanel(panel, -1, panel.drawingPanel)
        panel.fileControlPanel = FileControlPanel(panel, -1)
        panel.plotControlPanel = PlotControlPanel(panel, -1)

        # add control panels to vertical sizer
        controlbox = wx.FlexGridSizer(wx.VERTICAL)
        controlbox.Add(panel.movieControlPanel, 1, wx.EXPAND)
        controlbox.Add(panel.fileControlPanel, 1, wx.EXPAND)
        controlbox.Add(panel.plotControlPanel, 1, wx.EXPAND)
        controlbox.Add(panel.colorControlPanel, 1, wx.EXPAND)
        controlbox.AddGrowableRow(3)

        hbox = wx.BoxSizer(wx.HORIZONTAL)
        hbox.Add(panel.drawingPanel, 1, wx.EXPAND)
        hbox.Add(controlbox, 0, wx.EXPAND)

        panel.SetSizer(hbox)
        panel.SetAutoLayout(1)
        hbox.Fit(panel)

        self.Show(True)

    def OnClickRedraw(self, event):
        print "redraw at frame"
        self.drawingPanel.Refresh()


def main():
    app = wx.App(False)
    frame = MainFrame(None, -1, 'pygs: a python grana stack visualizer')
    app.MainLoop()
 
main()
