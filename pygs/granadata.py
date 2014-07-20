#!/usr/bin/env python

""" @module granadata
Contains classes for GranaData, Datapoint, Particle, and Params.
For use with pygs visualizer or as API
"""

import os
from math import *
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.delaunay as delaunay
import re
import csv
import fnmatch
from matplotlib import colors
from matplotlib import cm


class GranaData:
    """Class that stores, manages, and analyzes configuration data
    for multilayer 2D systems.
    """
    """
    @var trajectory Stores Datapoint objects in dict, usage trajectory[time][type][layer][particle]
    @var particles Stores Particle objects in list, usage particles[type][layer][particle]
    @var cluster_sizes Stores energy transfer clusters in dict, usage cluster_sizes[time][layer][clusterid][type]
    @var array_sizes Stores PSII arrays in dict, usage array_sizes[time][layer][clusterid][type]
    @var params Params object
    """

    def __init__(self):
        """public: constructor
        @param self The object pointer
        """
        self.trajectory = dict()
        self.particles = []
        for i_type in range(9):
            self.particles.append([]) # particle type list
            self.particles[i_type].append([]); # layer 0
            self.particles[i_type].append([]); #layer 1
        self.cluster_sizes = dict()
        self.array_sizes = dict()
        self.params = Params()

    def AddParticle(self, particle):
        """private: adds a Particle object to self.particles
        @param self The object pointer
        @param particle The Particle object to add
        """
        (type, layer, id) = particle.getID();
        particle_set = self.particles[type][layer]
        while len(particle_set) <= id:
            particle_set.append(False)
        particle_set[id] = particle

    def hasParticle(self, type, layer, id):
        """private: returns a bool indicating if a particle exists with the given identifiers
        @param self The object pointer
        @param type Int for particle type, eg 0 for LHCII or 1 for PSII
        @param layer Int for which layer particle is in, eg 0 or 1
        @param id Int for particle id number, ie index in list
        @retval bool True if present, false if not
        """
        if len(self.particles[type][layer]) <= id:
            return False
        else:
            if self.particles[type][layer][id]:
                return True
            else:
                return False

    def SetUpTimepoint(self, time, Lx=300, Ly=300, grana_rad=0, stroma_width=0):
        """private: sets up a new timepoint in self.trajectory, self.params, and self.cluster_sizes
        @param self The object pointer
        @param time Int for the new time
        @param L Box side length in x and y, default 300
        @param grana_rad Radius of grana disc, default 0
        @param stroma_width Width of stroma lamellae, default 0
        """
        self.trajectory[time] = []
        for i_type in range(len(self.particles)):
            self.trajectory[time].append([])
            type_list = self.particles[i_type]
            for i_layer in range(len(type_list)):
                self.trajectory[time][i_type].append([])
                layer = type_list[i_layer]
                for i_particle in layer:
                    self.trajectory[time][i_type][i_layer].append(False)
        self.cluster_sizes[time] = dict()
        self.params.addSystemSize(time, Lx, Ly, grana_rad, stroma_width)
        self.params.addTime(time)

    def AddPositionToTraj(self, time, type, layer, id, x, y, theta, energy, region=0, cluster=-1):
        """private: adds single particle/single time properties to self.trajectory

        by adding a new Datapoint object and filling missing indices with False
        @param self The object pointer
        @param time Int for time
        @param type Int for particle type, eg 0 for LHCII or 1 for PSII
        @param layer Int for which layer the particle is in, eg 0 or 1
        @param id Int for particle id number, ie index in list
        @param x Float for x position in 2D layer
        @param y Float for y position in 2D layer
        @param theta Float for rotation angle in radians
        @param energy Float for particle energy
        @param region Int indicating region within thylakoid, default 0 (0=grana, 1=stroma, 2=disallowed)
        @param cluster Int indicating cluster particle has been assigned to, default -1
        @retval bool True if file read successfully, False if not
        """
        particle_list = self.trajectory[time][type][layer]
        while len(particle_list) <= id:
            particle_list.append(False)
        particle_list[id] = Datapoint(x, y, theta, energy, region, cluster)
        return True

    def readFile(self, filename):
        """public: load data from file

        file name format: config_step[time]*L[boxsize]*.txt \n
        or config_step[time]*Lx[width]*Ly[height]*.txt

        file format: \n
        col1 (required) particle type, eg 0 for LHCII or 1 for PSII \n
        col2 (required) particle id, int starting at 0 \n
        col3 (required) particle x coordinate \n
        col4 (required) particle y coordinate \n
        col5 (required) layer id, eg 0 for bottom or 1 for top \n
        col6 (required) particle radius (should always be positive) \n
        col7 (required) particle length (0 for discs, positive for discorectangles) \n
        col8 (required) particle rotation in radians \n
        col9 (required) particle interaction state (aka phosphorylation), 0=interacting, 1=not interacting \n
        col10 (required) particle energy in kT \n
        col11 (optional) region id, default 0 (0=grana, 1=stroma, 2=disallowed) \n
        col12 (optional) cluster id if particle has been assigned to a cluster during simulation \n
        @param self The object pointer
        @param filename Relative or absolute path to file
        @retval True
        """
        basename = os.path.basename(filename)

        # get time from filename
        mtime = re.search('config_step(\d+)', basename)
        time = 0
        if mtime:
            time = int(mtime.group(1))
        else:
            print "can't find the time in filename", filename, "read as initial config"
        if time % self.params.stride != 0:
            return

        # get system size from filename
        mL = re.search('_L(\d+.\d+)', basename)
        mLx = re.search('_Lx(\d+.\d+)', basename)
        mLy = re.search('_Ly(\d+.\d+)', basename)
        Lx, Ly = 300, 300
        if mL:
            L = float(mL.group(1))
            Lx = L
            Ly = L
        else:
            if mLx and mLy:
                Lx = float(mLx.group(1))
                Ly = float(mLy.group(1))
            else:
                print "can't find any system size in filename", filename, "setting size to 300x300"

        # set up timepoint 
        print "reading from", filename, "into time", time, "with system size", Lx, "x", Ly
        if time not in self.trajectory.keys():
             self.SetUpTimepoint(time, Lx, Ly)
        #else:
        #    self.SetUpTimepoint(time)
        
        # read from file
        file = open(filename)
        nread = 0
        nadded = 0
        count = 0
        for line in file:
            # get data
            linedata = line.split()

            # read as a normal output file
            if (len(linedata) == 10):
                type = int(linedata[0])
                id = int(linedata[1])
                x = float(linedata[2])
                y = float(linedata[3])
                layer = int(linedata[4])
                radius = float(linedata[5])
                length = float(linedata[6])
                theta = float(linedata[7])
                phos = bool(int(linedata[8]))
                energy = float(linedata[9])
                region = 0
                cluster = -1

            # read as an output file for grana+stroma configs
            elif (len(linedata) == 11):
                type = int(linedata[0])
                id = int(linedata[1])
                x = float(linedata[2])
                y = float(linedata[3])
                layer = int(linedata[4])
                radius = float(linedata[5])
                length = float(linedata[6])
                theta = float(linedata[7])
                phos = bool(int(linedata[8]))
                energy = float(linedata[9])
                region = int(linedata[10])
                cluster = -1

            # read as an output file for grana+stroma configs w clusters
            elif (len(linedata) == 12):
                type = int(linedata[0])
                id = int(linedata[1])
                x = float(linedata[2])
                y = float(linedata[3])
                layer = int(linedata[4])
                radius = float(linedata[5])
                length = float(linedata[6])
                theta = float(linedata[7])
                phos = bool(int(linedata[8]))
                energy = float(linedata[9])
                region = int(linedata[10])
                cluster = int(linedata[11])
                
            # read LHCII from old input file format
            elif (int(linedata[0]) == 0 and len(linedata) == 3 and time == 0):
                type = int(linedata[0])
                id = count
                x = float(linedata[1])
                y = float(linedata[2])
                layer = 0
                radius = 6.5/2.0
                length = 0
                theta = 0
                phos = False
                energy = 0
                region = 0
                cluster = -1
                count = count + 1

            # read PSII from old input file format
            elif (int(linedata[0]) == 1 and len(linedata) == 4 and time == 0):
                type = int(linedata[0])
                id = count
                x = float(linedata[1])
                y = float(linedata[2])
                layer = 0
                radius = 6
                length = 26.5 - 12
                theta = float(linedata[3])
                phos = False
                energy = 0
                region = 0
                cluster = -1
                count = count + 1 

            # read disc from old input file format
            elif (int(linedata[0]) == 5 and len(linedata) == 7 and time == 0):
                type = int(linedata[0])
                id = count
                x = float(linedata[1])
                y = float(linedata[2])
                layer = 0
                radius = float(linedata[3])
                length = 0
                theta = 0
                phos = False
                energy = 0
                region = 0
                cluster = -1
                count = count + 1 

            else:
                print "got bad line", line
                return False

            # initialize particle if needed
            if not self.hasParticle(type, layer, id):
            #    print "adding particle", type, layer, id
                particle = Particle(type, id, layer, radius, length, phos)
                self.AddParticle(particle)
                nadded = nadded + 1

            # add position to trajectory
            if self.AddPositionToTraj(time, type, layer, id, x, y, theta, energy, region, cluster):
                nread = nread + 1
            else:
                print "problem with line", line
            
        # finish
        print "read", nread, "lines and added", nadded, "particles; have particles per layer", [ [len(plist) for plist in llist] for llist in self.trajectory[time] ]
        self.params.setTime(time)
        return True

    def PBC_rsq(self, x1, y1, x2, y2, time):
        """public: return distance-squared between two points, using PBC at given time
        @param self The object pointer
        @param x1 x coord of first point
        @param y1 y coord of first point
        @param x1 x coord of second point
        @param y1 y coord of second point
        @retval rsq double dx**2 + dy**2
        """
        dx = self.PBC_diff(x1, x2, self.params.getWidth(time))
        dy = self.PBC_diff(y1, y2, self.params.getHeight(time))
        return dx*dx + dy*dy
 
    def PBC_diff(self, val1, val2, L):
        """public: finds distance in 1D between points using periodic boundary conditions

        uses self.params.doPBC
        @param self The object pointer
        @param val1 First point
        @param val2 Second point
        @param L Box side length
        @retval Float distance between points
        """
        # calculate signed PBC distance for val1 - val2
        halfL = L * 0.5

        diff = val1 - val2
        if not self.params.doPBC:
            return diff
        elif diff > halfL:
            return diff - L
        elif diff > -halfL:
            return diff
        else:
            return diff + L

    def squared_dist_point_lineseg(self, xp, yp, x1, y1, x2, y2, length, Lx, Ly):
        """public: finds distance between point and line segment within a periodic box
        @param self The object pointer
        @param xp x coordinate of point
        @param yp y coordinate of point
        @param x1 x coordinate of endpoint 1
        @param y1 y coordinate of endpoint 1
        @param x2 x coordinate of endpoint 2
        @param y2 y coordinate of endpoint 2
        @param length length of line segment, avoids recalculating it
        @param Lx box side length in x dimension
        @param Ly box side length in y dimension
        @retval Float distance
        """
        line_seg_x = self.PBC_diff(x2, x1, Lx)
        line_seg_y = self.PBC_diff(y2, y1, Ly)
        projection = (self.PBC_diff(xp, x1, Lx) * line_seg_x + self.PBC_diff(yp, y1, Ly) * line_seg_y) / length / length
        closest_x = 0
        closest_y = 0
        if projection < 0:
            closest_x = x1
            closest_y = y1
        elif projection > 1:
            closest_x = x2
            closest_y = y2
        else:
            closest_x = self.PBC_diff(x1 + projection * line_seg_x, 0, Lx)
            closest_y = self.PBC_diff(y1 + projection * line_seg_y, 0, Ly)

        dx = self.PBC_diff(xp, closest_x, Lx)
        dy = self.PBC_diff(yp, closest_y, Ly)
        return dx*dx+dy*dy

    def getDelaunayNbhrs(self, itype, layer, time, do_plot=False):
        """private: calculate the nearest neighbors of all particles for a given type, layer, and time using Delaunay triangulation

        stores values in Datapoint.nbhrs
        @param self The object pointer
        @param type Int for particle type
        @param layer Int for layer id
        @param time Int for time
        @param do_plot Bool for whether to display Delaunay triangulation with matplotlib, True = do draw, False = don't draw
        @retval True if nbhrs were read from file, False if nbhrs were calculated from scratch
        """
        # https://github.com/matplotlib/matplotlib/blob/master/lib/matplotlib/tri/triangulation.py
        # http://www.astro.rug.nl/efidad/matplotlib.delaunay.triangulate.html#Triangulation
        
        # first look for neighbor file
        dfiletarget = "delaunay_type"+str(itype)+"_step"+str(time)+"_layer"+str(layer)+".txt"
        dfilename = False
        for file in os.listdir(self.params.output_path):
            if fnmatch.fnmatch(file, dfiletarget):
                dfilename = file
        # read file if found
        if dfilename:
            counts = dict()
            print "reading Delaunay neighbor list from file", dfilename
            dfilereader = csv.reader(open(self.params.output_path+"/"+dfilename, "r"), delimiter=' ')
            for row in dfilereader:
                ip = int(row.pop(0))
                if self.trajectory[time][itype][layer][ip]:
                    self.trajectory[time][itype][layer][ip].nbhrs = [int(ipn) for ipn in row]
                    # histogram number of neighbors
                    if len(self.trajectory[time][itype][layer][ip].nbhrs) in counts:
                        counts[len(self.trajectory[time][itype][layer][ip].nbhrs)] += 1
                    else:
                        counts[len(self.trajectory[time][itype][layer][ip].nbhrs)] = 1
            print "nbhr histogram:", counts, sum(counts.values())
            return True

        # if no preexisting file, calculate from scratch
        print "calculating Delaunay neighbor list for type, layer, time, PBC =", itype, layer, time, self.params.doPBC

        # collect x and y coords into list
        xcoords = []
        ycoords = []
        indices = []
        for ip in range(len(self.trajectory[time][itype][layer])):
            if self.trajectory[time][itype][layer][ip]:
                pos = self.trajectory[time][itype][layer][ip]
                xcoords.append(pos.x)
                ycoords.append(pos.y)
                indices.append(ip)
        npts = len(indices)

        # prepare periodically replicated coords if needed
        # order of periodic images:
        # 1 2 3
        # 4 0 5
        # 6 7 8
        if self.params.doPBC:
            for yimage in [1,0,-1]:
                for ximage in [-1,0,1]:
                    if ximage == 0 and yimage == 0:
                        continue
                    else:
                        for ip in indices:
                            xcoords.append(xcoords[ip] + ximage * self.params.getWidth(time))
                            ycoords.append(ycoords[ip] + yimage * self.params.getHeight(time))

        # use matplotlib's Delaunay triangulation
        # https://github.com/matplotlib/matplotlib/blob/master/lib/matplotlib/tri/triangulation.py
        # centers is Ntri x 2 float array of circumcenter x,y coords (ie vertices in Voronoi)
        # edges is N? x 2 int array of indices defining triangle edges
        # triangles is Ntri x 3 int array of triangle vertex indices, ordered counterclockwise
        # tri_nhbrs is Ntri x 3 int array of triangle indices that share edges (null = -1)
        #   tri_nhbrs[i,j] is the triangle that is the neighbor to the edge from
        #   point index triangles[i,j] to point index triangles[i,(j+1)%3].
        centers, edges, triangles, tri_nbhrs = delaunay.delaunay(xcoords, ycoords)

        # if PBC, use own code
        #    nbhrs, triangles = self.delaunayPBC(xcoords, ycoords, time)
        
        # store triangles as set of tuples in triangulation
        triangulation = set()
        for tri in triangles:
            sorted_tri = tuple(sorted(tri))
            if self.params.doPBC and sorted_tri[1] >= npts:
                # if PBC, don't store tri if 2 or more vertices outside original image
                #print "rejecting tri", sorted_tri
                continue
            else:
                triangulation.add(sorted_tri)

        # identify nearest neighbors in triangulation
        # nbhrs is a dict with key=vertex, val=set of neighbor tmp particle ids
        # all particle indices are in original image
        nbhrs = dict(zip(indices, [set() for i in indices]))
        for tri in triangulation:
            imaged_tri = [ip % npts for ip in tri]
            for ivert in range(3): # loop through triangle vertices
                for offset in range(2): # loop through the two neighboring vertices
                    nbhrs[imaged_tri[ivert]].add(imaged_tri[(ivert+offset)%3])

        # add neighbors to particles in self.trajectory
        counts = dict()
        print "writing Delaunay neighbor list to file", dfiletarget
        dfilewriter = csv.writer(open(self.params.output_path+"/"+dfiletarget, "wb"), delimiter=' ')
        for real_ip in range(len(self.trajectory[time][itype][layer])):
            if self.trajectory[time][itype][layer][real_ip]:
                # store in nbhrs
                tmp_ip = indices.index(real_ip)
                self.trajectory[time][itype][layer][real_ip].nbhrs = [ indices[nbhr_ip] for nbhr_ip in nbhrs[tmp_ip] ]
                
                # write to file
                dfilewriter.writerow([real_ip] + self.trajectory[time][itype][layer][real_ip].nbhrs)
                
                # histogram number of neighbors
                if len(nbhrs[tmp_ip]) in counts:
                    counts[len(nbhrs[tmp_ip])] += 1
                else:
                    counts[len(nbhrs[tmp_ip])] = 1
                        
       # print sorted(list(triangulation))
       # print list(triangulation)
        print len(triangulation), len(triangles)
        print "nbhr histogram:", counts, sum(counts.values())
            
        # plot triangles
        if do_plot:
            npx = np.asarray(xcoords)
            npy = np.asarray(ycoords)
            for t in triangulation:
                t_i = [t[0], t[1], t[2], t[0]]
                x_i = list(npx[t_i])
                y_i = list(npy[t_i])
                
                if False: #self.params.doPBC:
                    min_xi = x_i.index(min(x_i)) #np.where(x_i == min(x_i))
                    max_xi = x_i.index(max(x_i)) #np.where(x_i == max(x_i))[0]
                    med_xi = (set([0,1,2]) - set([min_xi, max_xi])).pop()
                    if x_i[max_xi] - x_i[med_xi] != self.PBC_diff(x_i[max_xi], x_i[med_xi], self.params.getWidth(time)):
                        x_i[max_xi] -= self.params.getWidth(time)
                    elif x_i[med_xi] - x_i[min_xi] != self.PBC_diff(x_i[med_xi], x_i[min_xi], self.params.getWidth(time)):
                        x_i[min_xi] += self.params.getWidth(time) 
                    elif x_i[max_xi] - x_i[min_xi] != self.PBC_diff(x_i[max_xi], x_i[min_xi], self.params.getWidth(time)):
                        if self.params.getWidth(time) - x_i[max_xi] > x_i[min_xi]:
                            x_i[min_xi] += self.params.getWidth(time)
                        else:
                            x_i[max_xi] -= self.params.getWidth(time)

                    x_i[3] = x_i[0]
                
                    min_yi = y_i.index(min(y_i)) #np.where(npy == min(y_i))
                    max_yi = y_i.index(max(y_i)) #np.where(npy == max(y_i))
                    med_yi = (set([0,1,2]) - set([min_yi, max_yi])).pop()
                    if y_i[max_yi] - y_i[med_yi] != self.PBC_diff(y_i[max_yi], y_i[med_yi], self.params.getHeight(time)):
                        y_i[max_yi] -= self.params.getHeight(time)
                    elif y_i[med_yi] - y_i[min_yi] != self.PBC_diff(y_i[med_yi], y_i[min_yi], self.params.getHeight(time)):
                        y_i[min_yi] += self.params.getHeight(time) 
                    elif y_i[max_yi] - y_i[min_yi] != self.PBC_diff(y_i[max_yi], y_i[min_yi], self.params.getHeight(time)):
                        if self.params.getHeight(time) - y_i[max_yi] > y_i[min_yi]:
                            y_i[min_yi] += self.params.getHeight(time)
                        else:
                            y_i[max_yi] -= self.params.getHeight(time)

                    y_i[3] = y_i[0]
                        
                plt.plot(x_i, y_i)
            plt.show()

        return False

    def delaunayPBC(self, xcoords, ycoords, time):
        """private: DEPRECATED helper function for GranaData.getDelaunayNbhrs
        @param self The object pointer
        @param xcoords List with x positions of all particles
        @param ycoords List with y positions of all particles
        @param time Int for time (needed for PBC distance calculations)
        @retval nbhrs Dict with key=vertex, val=set of neighbor tmp particle ids
        """
        # DOI 10.1103/PhysRevE.80.041101
        # http://www.cs.berkeley.edu/~jrs/meshpapers/SuThesis.ps.gz Chapter 2 sec 2.5, fig 2.13

        # q_edges and finished_edges are sets of edges
        # edges are sorted tuples (ip1, ip2) of tmp particle ids (actual particle ids are indices[ip])
        q_edges = set()
        finished_edges = set()
        # triangulation is a set of (ip1, ip2, ip3) vertices
        triangulation = set()
        # nbhrs is a dict with key=vertex, val=set of neighbor tmp particle ids
        nbhrs = dict()
        
        # find valid Delaunay triangle that includes first point
        # start by finding nearest neighbor
        minrsq = 10000000
        ip1 = 0
        ip2 = ip1
        npts = len(xcoords)
        for trial_ip2 in range(npts):
            if ip1 != trial_ip2:
                rsq = self.PBC_rsq(xcoords[ip1], ycoords[ip1], xcoords[trial_ip2], ycoords[trial_ip2], time)
                if rsq < minrsq:
                    minrsq = rsq
                    ip2 = trial_ip2
        # add first triangle
        found_tri = False
        for ip3 in range(len(xcoords)):
            if ip3 in [ip1, ip2]:
                continue
            found_tri = self.isDelaunayTriangle(ip1, ip2, ip3, xcoords, ycoords, time)
            if found_tri:
                break
        q_edges.add(tuple(sorted([ip1,ip2])))
        q_edges.add(tuple(sorted([ip2,ip3])))
        q_edges.add(tuple(sorted([ip1,ip3])))
        triangulation.add(tuple(sorted([ip1, ip2, ip3])))
        nbhrs[ip1] = set([ip2, ip3])
        nbhrs[ip2] = set([ip1, ip3])
        nbhrs[ip3] = set([ip2, ip1])
        print triangulation

        # build rest of triangulation
        while len(q_edges) > 0:
            # get random unfinished edge
            base_edge = q_edges.pop()
            # mark edge as finished
            finished_edges.add(base_edge)
        #    print "processing edge", base_edge
            # find vertices of already-found triangle(s) by intersection of nbhr sets
            (ip1, ip2) = base_edge
            shared_verts = nbhrs[ip1] & nbhrs[ip2]
  
            if len(shared_verts) > 0:

                # have at least one complete triangle, need to find other
                n_tris = 0

                # set up ordered list of possible third vertices
                close_verts = nbhrs[ip1] | nbhrs[ip2]
                for ip in shared_verts:
                    close_verts |= nbhrs[ip]
                all_verts = set(range(len(xcoords)))
                far_verts = all_verts - close_verts - set([ip1, ip2])
                verts_to_check = list(close_verts- set([ip1, ip2])) + list(far_verts)

                for ip3 in verts_to_check:

                    test_tri = tuple(sorted([ip1, ip2, ip3]))

                    # check if triangle has been found
                    if test_tri in triangulation:
                        # if triangle has been found before, move on
                        n_tris += 1

                    else:
                        # if triangle not yet found, check if it's valid
                        noip3_verts_to_check = verts_to_check[:]
                        noip3_verts_to_check.remove(ip3)
                        is_tri = self.isDelaunayTriangle(ip1, ip2, ip3, xcoords, ycoords, time, noip3_verts_to_check)

                        if is_tri:
                            # add new triangle
                            n_tris += 1
                            triangulation.add(test_tri)

                            # add new edges and nbhrs
                            if ip3 in nbhrs:
                                nbhrs[ip3].update([ip1, ip2])
                            else:
                                nbhrs[ip3] = set([ip1, ip2])
                            for vertex in base_edge:
                                new_edge = tuple(sorted([vertex,ip3]))
                                nbhrs[vertex].add(ip3)
                                if new_edge not in q_edges | finished_edges:
                                    q_edges.add(new_edge)
                                    if new_edge == (1, 16):
                                        print "adding edge", new_edge, base_edge, shared_verts

                        else:
                            # if not triangle, move on
                            continue

                    if n_tris == 2:
                        # found both triangles that share this edge, done
                        break

                # end check for triangles
                if self.params.doPBC and n_tris != 2:
                    print "WARNING: didn't find two triangles for edge", base_edge, shared_verts, n_tris, len(verts_to_check)
                    for ip in list(base_edge) + list(shared_verts):
                        print ip, xcoords[ip], ycoords[ip]
            else:
                print "WARNING: no nbhrs for edge", base_edge, finished_verts

        # check for valid triangulation (PBC only)
        if self.params.doPBC and len(triangulation) != 2*npts:
            print "WARNING: bad triangulation", len(triangulation), 2*npts, triangulation
            tri_copy = triangulation.copy()
            for edge in finished_edges:
                for third in nbhrs[edge[0]] & nbhrs[edge[1]]:
                    if tuple(sorted([edge[0], edge[1], third])) not in triangulation:
                        print "missing triangle is", edge, tuple(sorted([edge[0], edge[1], third]))
            raise

        # if good, return triangulation
        else:
            return nbhrs, triangulation


    def isDelaunayTriangle(self, ip1, ip2, ip3, xcoords, ycoords, time, verts_to_check = []):
        """private: helper function for GranaData.delaunayPBC
        @param self The object pointer
        @param ip1 Int for id of vertex 1 of facet
        @param ip2 Int for id of vertex 2 of facet
        @param ip3 Int for id of possible vertex 3
        @param xcoords List with x positions of all particles
        @param ycoords List with y positions of all particles
        @param time Int for time (needed for PBC distance calculations)
        @param nearby_ids (optional) set of ids to check first
        @retval bool True if ip3 completes a Delaunay triangle, False if not
        """

        x2 = self.PBC_diff(xcoords[ip2], xcoords[ip1], self.params.getWidth(time))
        y2 = self.PBC_diff(ycoords[ip2], ycoords[ip1], self.params.getHeight(time))
        
        # set up list of vertices to check, if not given
        if len(verts_to_check) == 0:
            all_verts = set(range(len(xcoords)))
            verts_to_check = list(all_verts- set([ip1, ip2, ip3]))

        # transform ip3 to coords where (x1,y1) = (0,0)
        x3 = self.PBC_diff(xcoords[ip3], xcoords[ip1], self.params.getWidth(time))
        y3 = self.PBC_diff(ycoords[ip3], ycoords[ip1], self.params.getHeight(time))

        # arrange CCW with Graham scan
        if x2*y3 - y2*x3 > 0:
            ccw_x2, ccw_y2 = x2, y2
            ccw_x3, ccw_y3 = x3, y3
            if set([149, 214, 207]) == set([ip1, ip2, ip3]):
                print "ccw2 = x2 at", ccw_x2, ccw_y2, ", ccw3 = x3 at", ccw_x3, ccw_y3
        else:
            ccw_x2, ccw_y2 = x3, y3
            ccw_x3, ccw_y3 = x2, y2
            if set([149, 214, 207]) == set([ip1, ip2, ip3]):
                print "ccw2 = x3, ccw3 = x2", len(verts_to_check)

        # check if any other point is inside circumcircle
        for ip4 in verts_to_check:
            # transform ip4 to coords where (x1,y1) = (0,0)
            x4 = self.PBC_diff(xcoords[ip4], xcoords[ip1], self.params.getWidth(time))
            y4 = self.PBC_diff(ycoords[ip4], ycoords[ip1], self.params.getHeight(time))
            
            # compute determinant
            mat = np.array([
                    [0, 0, 0, 1],
                    [ccw_x2, ccw_y2, ccw_x2*ccw_x2 + ccw_y2*ccw_y2, 1],
                    [ccw_x3, ccw_y3, ccw_x3*ccw_x3 + ccw_y3*ccw_y3, 1],
                    [x4, y4, x4*x4 + y4*y4, 1],
                    ])
            det = np.linalg.det(mat)

            if det > 0: # not a Delaunay triangle
                if set([149, 214, 207]) == set([ip1, ip2, ip3]):
                    print ip1, ip2, ip3, ip4, x4, y4, det
                return False
            
        # if we got outside the loop, ip3 is a valid neighbor
        return True

    def calculateDelaunayShapes(self, itype, id, layer, time):
        """private: calculate the local shape-based order parameter of a given particle type

        returns value and stores it in Particle.
        @param self The object pointer
        @param type Int for particle type
        @param id Int for particle id
        @param layer Int for layer id
        @param time Int for time
        @retval Float value of the hexatic order parameter, between 0 and 1
        """
        # use Lester's trig-free method for hexatic order parameter sum(exp(6 * i * theta))
       #     realpart = realpart + pow(c, 6) - 15.0 * pow(c, 4) * pow(s, 2) + 15.0 * pow (c, 2) * pow(s, 4) - pow(s, 6)
       #     impart = impart + 6.0 * pow(c, 5) * s - 20.0 * pow(c, 3) * pow(s, 3) + 6.0 * c * pow(s, 5)
       
        # better, use Aaron's shape matching method: http://www.annualreviews.org/doi/full/10.1146/annurev-conmatphys-062910-140526
        p1 = self.particles[itype][layer][id]
        pos1 = self.trajectory[time][itype][layer][id]

        if not pos1:
            return 0
        if pos1.hexatic:
            return pos1.hexatic

        typelist = [itype]
        if type in [0, 2]:
            typelist = [0, 2]

        # find neighbors if needed
        if len(pos1.nbhrs) == 0:
            self.getDelaunayNbhrs(itype, layer, time)
        
        # set up ideal crystals with shape, vector magnitude, weights, and cutoff value
        ideal_shape = dict()
        mags = dict()
        cutoff = dict()
        weights = dict()

        ideal_shape["c2s2boekema"] = dict(zip([7, 9, 10, 12], [0.12, 0.12, 0.21, 0.84]))
        mags["c2s2boekema"] = sqrt(sum([val*val for val in ideal_shape["c2s2boekema"].values()]))
        cutoff["c2s2boekema"] = 0.80

        ideal_shape["c2s2mboekema"] = dict(zip([4, 5, 6, 7, 8, 9, 10], [0.35, 0.05, 0.30, 0.05, 0.12, 0.10, 0.95]))
        mags["c2s2mboekema"] = sqrt(sum([val*val for val in ideal_shape["c2s2mboekema"].values()]))
        cutoff["c2s2mboekema"] = 0.75

        ideal_shape["c2s2m2kourilA"] = dict(zip([2, 3, 4, 5, 6, 7, 8, 9, 10], [0.15, 0.02, 0.38, 0.02, 0.75, 0.05, 0.32, 0.07, 0.87]))
        mags["c2s2m2kourilA"] = sqrt(sum([val*val for val in ideal_shape["c2s2m2kourilA"].values()]))
        cutoff["c2s2m2kourilA"] = 0.65

        ideal_shape["c2s2m2kourilB"] = dict(zip([2, 3, 4, 5, 6, 7, 8, 9, 10], [0.21, 0.0, 0.45, 0.02, 0.60, 0.03, 0.55, 0.01, 0.75]))     
        mags["c2s2m2kourilB"] = sqrt(sum([val*val for val in ideal_shape["c2s2m2kourilB"].values()]))
        cutoff["c2s2m2kourilB"] = 0.65
       
        # set up containers
        ls = ideal_shape[self.params.crystal].keys()
        zeros = [0 for l in ls]
        psis = dict(zip(ls, zeros))
        cosines = dict(zip(ls, zeros))
        sines = dict(zip(ls, zeros))
        nbhrs = len(pos1.nbhrs)        

        # calculate local shape contribution for each neighbor
        for ip in pos1.nbhrs:
            p2 = self.particles[type][layer][ip]
            if p2:
                pos2 = self.trajectory[time][type][layer][p2._id]
                dx = self.PBC_diff(pos1.x, pos2.x, self.params.getWidth(time))
                dy = self.PBC_diff(pos1.y, pos2.y, self.params.getHeight(time))
                r = sqrt(dx*dx + dy*dy)
                c = dx / r
                s = dy / r
                theta = atan2(s,c)
                for l in ls:
                    cosines[l] += cos(l*theta)
                    sines[l] += sin(l*theta)

        # calculate local shape psi(l)
        for l in ls:
            if nbhrs > 0:
                psis[l] = sqrt((cosines[l])**2 + (sines[l])**2) / nbhrs
            else:
                psis[l] = 0.0
        psimag = sqrt(sum([val*val for val in psis.values()]))
       
        # compare to reference crystal
     #   diff_shape = [abs(ideal_shape[self.params.crystal][l] - psi[l]) * weights[self.params.crystal][l] for l in ideal_shape[self.params.crystal].keys()]
        dot_shape = [ideal_shape[self.params.crystal][l] * psis[l] / mags[self.params.crystal] / psimag for l in ideal_shape[self.params.crystal].keys()]
        if sum(dot_shape) > cutoff[self.params.crystal]:
            pos1.hexatic = 0.95
        else:
            pos1.hexatic = 0.2
     #   pos1.hexatic = sum(dot_shape)
     #   pos1.hexatic = psis[4]
     #   if id == 20:
     #       print pos1.hexatic
     #   print id, pos1.hexatic
        return pos1.hexatic

    def calculateNematic(self, type, id, layer, time, returnnbhrs = False):
        """deprecated: calculate the local nematic order parameter of a given particle

        Use buildArrayClusters instead.

        only counts neighbors whose centers are within hard-coded distance of the given particle's center, ie 29.15 nm = 1.1*(length+2*radius)
        @param self The object pointer
        @param type Int for particle type
        @param id Int for particle id
        @param layer Int for layer id
        @param time Int for time
        @param returnnbhrs Bool for whether to return number of neighbors found, default False
        @retval Float value of the nematic order parameter (between 0 and 1), or number of neighbors
        """
        # local nematic order parameter: average 1/2 (3 cos^2(theta) - 1), where theta = angle diff btw neighbors
        p1 = self.particles[type][layer][id]
        pos1 = self.trajectory[time][type][layer][id]

        if not pos1:
            return 0
        if pos1.nematic:
            if returnnbhrs:
                return pos1.nbhrs
            else:
                return pos1.nematic

        typelist = [type]
        if type in [0, 2]:
            typelist = [0, 2]

        r_min = 2.0*p1._radius*1.5#0 # 6*2*1.3 = 15.6 nm
        r_max = (p1._length + 2.0*p1._radius) * 1.1 # (14.5+12)*1.1 = 29.15 nm
        nbhrs = 0.0
        sum = 0.0
        nbhrs_above_cutoff = 0.0
        close_nbhrs = 0.0
        for itype in typelist:
            for p2 in self.particles[itype][layer]:
                if p2:
                    if p2._id != id:
                        pos2 = self.trajectory[time][itype][layer][p2._id]
                        dx = self.PBC_diff(pos1.x, pos2.x, self.params.getWidth(time))
                        dy = self.PBC_diff(pos1.y, pos2.y, self.params.getHeight(time))
                        r = sqrt(dx*dx + dy*dy)
                        if r < r_max:# and r > r_min:
                            nbhrs = nbhrs + 1.0
                            theta = pos1.theta - pos2.theta
                            costheta = cos(theta)
                            val = (3.0*costheta*costheta - 1.0)/2.0
                            sum += val
                            if val > 0.9:
                                nbhrs_above_cutoff += 1.0
                                if r < r_min:
                                    close_nbhrs += 1.0

        if nbhrs == 0:
            nematic = 0
        else:
            nematic = sum / nbhrs
        self.trajectory[time][type][layer][id].nematic = nematic
        self.trajectory[time][type][layer][id].nbhrs = nbhrs
        if returnnbhrs:
            return (nbhrs_above_cutoff/6.0)*(close_nbhrs/2.0)
        else:
            return nematic

    def calculateHexatic(self, type, id, layer, time):
        # use Lester's trig-free method for hexatic order parameter sum(exp(6 * i * theta))                                  
        p1 = self.particles[type][layer][id]
        pos1 = self.trajectory[time][type][layer][id]

        if not pos1:
            return 0
      #  if hasattr(pos1, "hexatic"):
       #     return pos1.hexatic

        typelist = [type]
        if type in [0, 2]:
            typelist = [0, 2]

        r_cut = p1._radius * 2.0 * 1.73 # ~ cos(30/180*pi) * 2
        nbhrs = 0.0
        realpart = 0.0
        impart = 0.0
        for itype in typelist:
            for p2 in self.particles[itype][layer]:
                if p2:
                    if p2._id != id:
                        pos2 = self.trajectory[time][itype][layer][p2._id]
                        dx = self.PBC_diff(pos1.x, pos2.x, self.params.getWidth(time))
                        dy = self.PBC_diff(pos1.y, pos2.y, self.params.getHeight(time))
                        r = sqrt(dx*dx + dy*dy)
                        if r < r_cut:
                            nbhrs = nbhrs + 1.0
                            c = dx / r
                            s = dy / r
                            realpart = realpart + pow(c, 6) - 15.0 * pow(c, 4) * pow(s, 2) + 15.0 * pow (c, 2) * pow(s, 4) - pow(s, 6)
                            impart = impart + 6.0 * pow(c, 5) * s - 20.0 * pow(c, 3) * pow(s, 3) + 6.0 * c * pow(s, 5)

        if nbhrs < 2:
            hexatic = 0
        else:
            hexatic = 1.0 / nbhrs / nbhrs * (realpart*realpart + impart*impart)
        pos1.hexatic = hexatic
        return hexatic

    def buildETClusters(self, time, ilayer, rsqmax=2.0*2.0):
        """public: builds clusters using energy transfer criteria, stores info in self.cluster_sizes
        @param self The object pointer
        @param time Int for time
        @param ilayer Int for layer id
        @param rsqmax Float for maximum distance to be counted as energy transfer partners, default 4.0 nm
        """
        # based on Carl's code
        icluster = 0
        nbhrs_to_check = [] # Carl's "reserve"
        if time not in self.cluster_sizes:
            self.cluster_sizes[time] = dict()
        if ilayer not in self.cluster_sizes[time]:
            self.cluster_sizes[time][ilayer] = []

        maybetypes = [0, 1]
        types = []
        for itype in range(len(maybetypes)):
            if len(self.trajectory[time][itype]) > 0:
                types.append(maybetypes[itype])
        unsorted_clusters =[]

        # remove any remembered clusters
        #for itype in types:
        #    for ipos in range(len(self.trajectory[time][itype][ilayer])):
        #        self.trajectory[time][itype][ilayer][ipos].already_clustered = False

        # build cluster list
        for itype in types:
            for ipos in range(len(self.trajectory[time][itype][ilayer])):
                if self.trajectory[time][itype][ilayer][ipos]:
                    if not self.trajectory[time][itype][ilayer][ipos].already_clustered:
                        self.trajectory[time][itype][ilayer][ipos].already_clustered = True
                        nbhrs_to_check.append([itype, ipos])
                        unsorted_clusters.append([0 for i in range(max(types)+1)])
                    
                        # while there are nbhrs whose nbhrs haven't been checked yet...
                        while len(nbhrs_to_check) > 0:
                            # add reference particle to cluster
                            [itype1, ipos1] = nbhrs_to_check.pop()
                            self.trajectory[time][itype1][ilayer][ipos1].clusterID = icluster
                            unsorted_clusters[icluster][itype1] += 1
                            
                            # loop to find neighbors
                            for itype2 in types:
                                for ipos2 in range(len(self.trajectory[time][itype2][ilayer])):
                                    if self.trajectory[time][itype2][ilayer][ipos2]:
                                        if not self.trajectory[time][itype2][ilayer][ipos2].already_clustered:
                                            sep = self.getClosestApproach(time, ilayer, itype1, ipos1, itype2, ipos2)
                                            if sep*sep < rsqmax:
                                                # add actual neighbor to queue
                                                nbhrs_to_check.append([itype2, ipos2])
                                                self.trajectory[time][itype2][ilayer][ipos2].already_clustered = True
                        # finished this cluster, on to the next one
                        icluster += 1

        # sort cluster list and re-id
        sortedIDs = sorted([[val, i] for i, val in enumerate(unsorted_clusters)], reverse=True)
        IDmap = []
        for i, val in enumerate(unsorted_clusters):
            newID = sortedIDs.index([val, i])
            IDmap.append(newID)
        for itype in types:
            for ipos in range(len(self.trajectory[time][itype][ilayer])):
                if self.trajectory[time][itype][ilayer][ipos]:
                    if self.trajectory[time][itype][ilayer][ipos].already_clustered:
                        self.trajectory[time][itype][ilayer][ipos].clusterID = IDmap[self.trajectory[time][itype][ilayer][ipos].clusterID]
        self.cluster_sizes[time][ilayer] = [val for [val, i] in sortedIDs]
        print "found", len(self.cluster_sizes[time][ilayer]), "clusters at time", time, "in layer", ilayer#, "with sizes:", self.cluster_sizes[time][ilayer]

    def buildArrayClusters(self, time, ilayer):
        """public: builds clusters using PSII semi-crystalline array criteria, stores info in self.array_sizes
        @param self The object pointer
        @param time Int for time
        @param ilayer Int for layer id
        @retval largest_cluster Int for size of largest cluster
        """
        # based on Carl's code
        iarray = 0
        nbhrs_to_check = [] # Carl's "reserve"
        far_nbhrs_to_check = [] # secondary queue
        if time not in self.array_sizes:
            self.array_sizes[time] = dict()
        if ilayer not in self.array_sizes[time]:
            self.array_sizes[time][ilayer] = []

        types = [1]
        unsorted_arrays =[]

        r_min = 12.0*1.5#0 # 6*2*1.3 = 15.6 nm
        r_max = 26.5 * 1.1 # (14.5+12)*1.1 = 29.15 nm

        perp_cutoff = 14.0
        par_cutoff = 26.5

        # remove any remembered arrays
        #for itype in types:
        #    for ipos in range(len(self.trajectory[time][itype][ilayer])):
        #        self.trajectory[time][itype][ilayer][ipos].already_arrayed = False

        # build array list
        for itype in types:
            for ipos in range(len(self.trajectory[time][itype][ilayer])):
                #print "checking particle", ipos
                if self.trajectory[time][itype][ilayer][ipos]:
                    if not self.trajectory[time][itype][ilayer][ipos].already_arrayed:
                        #print ipos, "not yet arrayed so checking for nbhrs"
                        self.trajectory[time][itype][ilayer][ipos].already_arrayed = True
                        nbhrs_to_check.append([itype, ipos])
                        unsorted_arrays.append([0 for i in range(max(types)+1)])
                    
                        # while there are nbhrs whose nbhrs haven't been checked yet...
                        while len(nbhrs_to_check) > 0:
                            # add reference particle to array
                            [itype1, ipos1] = nbhrs_to_check.pop()
                            #print "checking particle", ipos1, "for nbhrs, started with", ipos
                            #if self.trajectory[time][itype1][ilayer][ipos1].already_arrayed:
                                #print "oops", ipos1, "is already in array", self.trajectory[time][itype1][ilayer][ipos1].arrayID
                            self.trajectory[time][itype1][ilayer][ipos1].arrayID = iarray
                            unsorted_arrays[iarray][itype1] += 1
                            #print iarray, unsorted_arrays[iarray]

                            # set up secondary nbhrs_to_check list, for ipos1 particle that isn't yet in linear array
                            far_nbhrs_to_check = []
                            
                            # loop to find neighbors
                            for itype2 in types:
                                for ipos2 in range(len(self.trajectory[time][itype2][ilayer])):
                                    if self.trajectory[time][itype2][ilayer][ipos2]:
                                        if not self.trajectory[time][itype2][ilayer][ipos2].already_arrayed:
                                            #print "checking particle", ipos2, "for nbhr to", ipos1

                                            # get pairwise nematic value
                                            pos1 = self.trajectory[time][itype1][ilayer][ipos1]
                                            pos2 = self.trajectory[time][itype2][ilayer][ipos2]
                                            theta = pos1.theta - pos2.theta
                                            costheta = cos(theta)
                                            val = (3.0*costheta*costheta - 1.0)/2.0

                                            # if pairwise value above cutoff, check distance
                                            if val > 0.9:
                                                # get distance between centers
                                                dx = self.PBC_diff(pos1.x, pos2.x, self.params.getWidth(time))
                                                dy = self.PBC_diff(pos1.y, pos2.y, self.params.getHeight(time))

                                                if True: # new method
                                                    unit_vector = [cos(pos1.theta), sin(pos1.theta)]
                                                    proj = unit_vector[0]*dx + unit_vector[1]*dy
                                                    par_dist = abs(proj)
                                                    perp_dist = sqrt(dx*dx + dy*dy - par_dist*par_dist)
                                                   # print "proj", proj, "par", par_dist, "perp", perp_dist, "unit", unit_vector
                                                    if par_dist < par_cutoff:
                                                        if perp_dist < perp_cutoff:
                                                            #print "in cutoffs, adding", ipos2
                                                            # add neighbor to queue
                                                            nbhrs_to_check.append([itype2, ipos2])
                                                            self.trajectory[time][itype2][ilayer][ipos2].already_arrayed = True
                                                            self.trajectory[time][itype2][ilayer][ipos2].in_linear_array = True
                                                else: # old method
                                                    r = sqrt(dx*dx + dy*dy)
                                                    if r < r_min:
                                                        # add close neighbor to close queue
                                                        nbhrs_to_check.append([itype2, ipos2])
                                                        self.trajectory[time][itype2][ilayer][ipos2].already_arrayed = True
                                                        self.trajectory[time][itype2][ilayer][ipos2].in_linear_array = True
                                                        
                                                        # if this is ipos1's first close neighbor, recheck its far nbhrs
                                                        if not self.trajectory[time][itype1][ilayer][ipos1].in_linear_array:
                                                            self.trajectory[time][itype1][ilayer][ipos1].in_linear_array = True
                                                            nbhrs_to_check.extend(far_nbhrs_to_check)
                                                            far_nbhrs_to_check = []
                                                    elif r < r_max:
                                                        if self.trajectory[time][itype1][ilayer][ipos1].in_linear_array:
                                                            # add far neighbor to queue, if ipos1 particle is in linear array
                                                            nbhrs_to_check.append([itype2, ipos2])
                                                            self.trajectory[time][itype2][ilayer][ipos2].already_arrayed = True
                                                        else:
                                                            # add far nbhr to tmp queue
                                                            far_nbhrs_to_check.append([itype2, ipos2])

                        # finished this array, on to the next one
                        iarray += 1

        # sort array list and re-id
        sortedIDs = sorted([[val, i] for i, val in enumerate(unsorted_arrays)], reverse=True)
        IDmap = []
        for i, val in enumerate(unsorted_arrays):
            newID = sortedIDs.index([val, i])
            IDmap.append(newID)
        for itype in types:
            for ipos in range(len(self.trajectory[time][itype][ilayer])):
                if self.trajectory[time][itype][ilayer][ipos]:
                    if self.trajectory[time][itype][ilayer][ipos].already_arrayed:
                        self.trajectory[time][itype][ilayer][ipos].arrayID = IDmap[self.trajectory[time][itype][ilayer][ipos].arrayID]
        self.array_sizes[time][ilayer] = [val for [val, i] in sortedIDs]
        total_ps = 0
        for arr in self.array_sizes[time][ilayer]:
            total_ps += arr[1]
        print "found", len(self.array_sizes[time][ilayer]), "arrays at time", time, "in layer", ilayer, "with largest cluster:", self.array_sizes[time][ilayer][0], "and total PSII", total_ps
        return self.array_sizes[time][ilayer][0]

    def getClosestApproach(self, time, ilayer, itype1, ipos1, itype2, ipos2):
        """private: returns closest distance between *edges* of particles, not centers
        @param self The object pointer
        @param time Int for time
        @param ilayer Int for layer id
        @param itype1 Int for type of particle 1
        @param ipos1 Int for particle id of particle 1
        @param itype2 Int for type of particle 2
        @param ipos2 Int for particle id of particle 2
        @retval Float distance
        """
        # set up box size
        Lx = self.params.getWidth(time)
        Ly = self.params.getHeight(time)

        # set up first particle
        pos1 = self.trajectory[time][itype1][ilayer][ipos1]
        p1 = self.particles[itype1][ilayer][ipos1]
            
        # set up second particle
        pos2 = self.trajectory[time][itype2][ilayer][ipos2]
        p2 = self.particles[itype2][ilayer][ipos2]

        # do centers first
        dx = self.PBC_diff(pos1.x, pos2.x, Lx)
        dy = self.PBC_diff(pos1.y, pos2.y, Ly)
        min_rsq = dx*dx + dy*dy
        if p1._length > 0:
            self.trajectory[time][itype1][ilayer][ipos1].findEnds(p1._length / 2.0)
            
            if p2._length > 0:
                self.trajectory[time][itype2][ilayer][ipos2].findEnds(p2._length / 2.0)
                
                # distance between two rods
                rsq_p1end1_p2 = self.squared_dist_point_lineseg(pos1.end1_x, pos1.end1_y, pos2.end1_x, pos2.end1_y, pos2.end2_x, pos2.end2_y, p2._length, Lx, Ly)
                if rsq_p1end1_p2 < min_rsq: min_rsq = rsq_p1end1_p2
                rsq_p1end2_p2 = self.squared_dist_point_lineseg(pos1.end2_x, pos1.end2_y, pos2.end1_x, pos2.end1_y, pos2.end2_x, pos2.end2_y, p2._length, Lx, Ly)
                if rsq_p1end2_p2 < min_rsq: min_rsq = rsq_p1end2_p2
                rsq_p2end1_p1 = self.squared_dist_point_lineseg(pos2.end1_x, pos2.end1_y, pos1.end1_x, pos1.end1_y, pos1.end2_x, pos1.end2_y, p1._length, Lx, Ly)
                if rsq_p2end1_p1 < min_rsq: min_rsq = rsq_p2end1_p1
                rsq_p2end2_p1 = self.squared_dist_point_lineseg(pos2.end2_x, pos2.end2_y, pos1.end1_x, pos1.end1_y, pos1.end2_x, pos1.end2_y, p1._length, Lx, Ly)
                if rsq_p2end2_p1 < min_rsq: min_rsq = rsq_p2end2_p1
            else:
                # distance between p1 rod and p2 disc
                rsq_p1_p2 = self.squared_dist_point_lineseg(pos2.x, pos2.y, pos1.end1_x, pos1.end1_y, pos1.end2_x, pos1.end2_y, p1._length, Lx, Ly)
                if rsq_p1_p2 < min_rsq: min_rsq = rsq_p1_p2
        else:
            if p2._length > 0:
                # distance between p1 disc and p2 rod
                rsq_p1_p2 = self.squared_dist_point_lineseg(pos1.x, pos1.y, pos2.end1_x, pos2.end1_y, pos2.end2_x, pos2.end2_y, p2._length, Lx, Ly)
                if rsq_p1_p2 < min_rsq: min_rsq = rsq_p1_p2
            #  no other cases, since already computed disc-disc distance
        return sqrt(min_rsq) - p1._radius - p2._radius

    def TagOnBoundary(self, time, ilayer, itype, cell_size):
        """tag particles that are near a boundary or region of missing data, based on given distance (assuming no PBC)
        @param self The object pointer
        @param time Int for time
        @param ilayer Int for layer id
        @param itype Int for type of particle
        @param cell_size Float for size of cell that, if empty, counts as missing data
        @retval particle density in valid region (particles/square nm)
        """
        # set up box and cell sizes
        Lx = self.params.getWidth(time)
        Ly = self.params.getHeight(time)
        print "going to tag on boundary with Lx=", Lx, "and Ly=", Ly, "with cell size", cell_size
        try:
            n_cells_x = int(Lx/cell_size)
            cell_size_x = Lx/float(n_cells_x)
            n_cells_y = int(Ly/cell_size)
            cell_size_y = Ly/float(n_cells_y)
        except ZeroDivisionError:
            print Lx, Ly, cell_size, n_cells_x #, n_cells_y
            raise
        cell_contents = [ [ [] for iy in range(n_cells_y) ] for ix in range(n_cells_x) ]
        is_bndy_cell = [ [ False for iy in range(n_cells_y) ] for ix in range(n_cells_x) ]

        # put particle ids in cells, and untag
        n_particles = 0
        for ip, pos in enumerate(self.trajectory[time][itype][ilayer]):
            if pos:
                cell_x = int(floor(pos.x/cell_size_x) % n_cells_x)
                cell_y = int(floor(pos.y/cell_size_y) % n_cells_y)
                try:
                    cell_contents[cell_x][cell_y].append(ip)
                except IndexError:
                    print cell_x, pos.x, cell_y, pos.y
                    raise
                n_particles += 1
                self.particles[itype][ilayer][ip]._tagged = False
      
        # identify and tag in boundary cells
        n_bndy_cells = 0
        n_bndy_particles = 0
        n_empty_cells = 0
        for ix, xlist in enumerate(cell_contents):
            for iy, plist in enumerate(xlist):
                # if cell is empty, tag self and neighboring cells (w/o PBC)
                if len(plist) == 0:
                    n_empty_cells += 1
                     # loop over self and neighboring cells (w/o PBC)
                    for inx in [ix-1, ix, ix+1]: # [ix-2, ix-1, ix, ix+1, ix+2]:
                        if inx >= 0 and inx < n_cells_x:
                            for iny in [iy-1, iy, iy+1]: # [iy-2, iy-1, iy, iy+1, iy+2]:
                                if iny >= 0 and iny < n_cells_y:
                                    if not is_bndy_cell[inx][iny]:
                                        # set cell and particles as bndy
                                        n_bndy_cells += 1
                                        is_bndy_cell[inx][iny] = True
                                        for ip in cell_contents[inx][iny]:
                                            self.particles[itype][ilayer][ip]._tagged = True
                                            n_bndy_particles += 1

                 # set edge cell as bndy, even if it has particles
                if (ix <= 1 or iy <= 1 or ix >= n_cells_x-2 or iy >= n_cells_y-2) and (not is_bndy_cell[ix][iy]):
                    n_bndy_cells += 1
                    is_bndy_cell[ix][iy] = True
                    for ip in plist:
                        self.particles[itype][ilayer][ip]._tagged = True
                        n_bndy_particles += 1
                   
        # compute particles per square micron in non-boundary cells
        nonbndy_area = float(n_cells_x * n_cells_y - n_bndy_cells) * cell_size_x * cell_size_y
        nonbndy_particles = float(n_particles - n_bndy_particles)
        if nonbndy_area > 0 and nonbndy_particles > 20:
            print "interior and total density per um^2:", nonbndy_particles/nonbndy_area*1000000.0, 1000000.0*n_bndy_particles/(float(n_bndy_cells - n_empty_cells) * cell_size_x * cell_size_y)
            print "number of nonboundary particles:", nonbndy_particles
            return nonbndy_particles/nonbndy_area
        else:
            print "too much boundary, try a different cell size?"
            return 0.0

    def TagMtrimers(self, time):
        """public: tag free LHCIIs (particle type 0) that are closest to an occupied in-plane M-trimer binding site (particle type m) at the given time, and untag all others
        @param self The object pointer
        @param time Int for time
        @retval n_any_mtrimers Int for number of tagged M-trimers (singly or doubly bound)
        @retval n_double_mtrimers Int for number of doubly-bound M-trimers (ie bridging between two m-sites)
        @retval n_unbound_trimers Int for number of unbound LHCII (neither S nor M)
        @retval n_occupied_msites Int for number of occupied M-binding sites
        """
        
        # set up time and types
        self.params.setTime(time)
        itypel = 0
        itypem = 4
        n_any_mtrimers = 0
        n_double_mtrimers = 0
        n_total_trimers = 0
        n_occupied_msites = 0

        # loop over layers
        for i_layer in range(len(self.trajectory[time][itypel])):
            # set up data
            lhc_pos = self.trajectory[time][itypel][i_layer]

            # start by untagging all free LHCIIs
            for ilhc in range(len(self.particles[itypel][i_layer])):
                if self.particles[itypel][i_layer][ilhc]:
                    self.particles[itypel][i_layer][ilhc]._tagged = False
                    n_total_trimers += 1

            # loop over occupied m-sites
            for msite in self.trajectory[time][itypem][i_layer]:
                if msite:
                    if msite.energy < 0:
                        n_occupied_msites += 1

                        # if occupied, look for neighbor within r_cut nm (plus round-off)
                        r_cut = 0.66*6.5 + 0.0001

                        # loop over LHCIIs
                        for ilhc in range(len(self.particles[itypel][i_layer])):
                            if self.particles[itypel][i_layer][ilhc] and self.trajectory[time][itypel][i_layer][ilhc]:
                                xdist = self.PBC_diff(self.trajectory[time][itypel][i_layer][ilhc].x, msite.x, self.params.getWidth(time))
                                ydist = self.PBC_diff(self.trajectory[time][itypel][i_layer][ilhc].y, msite.y, self.params.getHeight(time))
                                dist = sqrt(xdist*xdist + ydist*ydist)

                                # tag as neighbor if close enough
                                if dist < r_cut:
                                    if self.particles[itypel][i_layer][ilhc]._tagged:
                                        n_double_mtrimers += 1
                                    else:
                                        n_any_mtrimers += 1
                                        self.particles[itypel][i_layer][ilhc]._tagged = True
                                    break

        return (n_any_mtrimers, n_double_mtrimers, n_total_trimers-n_any_mtrimers, n_occupied_msites)

    def PlotMtrimers(self):
        """public: compute and print trajectory of M-LHCII statistics

        output file name: self.params.output_path+"/mtrimers_start"+str(times[0])+"_end"+str(maxtime)+"_stride"+str(stride)+".dat"

        output file format: \n
        col1 time \n
        col2 number of all M-LHCII \n
        col3 number of bridging (doubly-bound) M-LHCII \n
        col4 number of unbound (not M) free LHCII \n
        col5 number of m-sites occupied by S-LHCII \n
        col6 size of largest cluster \n
        @param self The object pointer
        """

        # set up times and stride
        mintime = self.params.time
        times = self.params.getTimes()
        stride = self.params.stride
        maxtime = times[-1]

        ofilename = self.params.output_path+"/mtrimers_start"+str(times[0])+"_end"+str(maxtime)+"_stride"+str(stride)+".dat"
        print "writing to", ofilename
        ofilewriter = csv.writer(open(ofilename, "wb"), delimiter=' ')

        for time in times:
            print "calculating M-trimer for time", time
            n_any_mtrimers, n_double_mtrimers, n_unbound_trimers, n_occupied_msites = self.TagMtrimers(time)

            if time not in self.array_sizes:
                self.buildArrayClusters(time, 1)

            ofilewriter.writerow([ time, n_any_mtrimers, n_double_mtrimers, n_unbound_trimers, n_occupied_msites-n_any_mtrimers-n_double_mtrimers, self.array_sizes[time][1][0][1] ])

    def PlotHexatic(self, i_type):
        """public: compute and print histogram of local hexatic order parameters of type-0 particles

        output file name: self.params.output_path+"/hexatic_type"+str(type)+"_start"+str(times[0])+"_end"+str(maxtime)+"_stride"+str(stride)+".dat"

        output file format: \n
        col1 bin for hexatic order parameter \n
        col2 number of particles in bin \n
        @param self The object pointer
        """
      #  type = int(0)

        # set up times and stride
        mintime = self.params.time
        times = self.params.getTimes()
        stride = self.params.stride
        
        # set up lists for single-particle hexatic and config-average hexatic
        single_hexatic = []
        config_hexatic = []
        fraction_ordered = []
        byparticle_hexatic = dict()
        for pos in self.particles[i_type][0] + self.particles[i_type][1]:
            if pos:
                byparticle_hexatic[pos._id] = []
        print 'going to look at '+str(len(byparticle_hexatic))+' particles'

        # loop over all start times
        for time in times:
            if time % stride == 0 and time >= mintime:
                print "time is", time
                # loop over layers
                for i_layer in range(len(self.trajectory[time][i_type])):
                    startlayer = self.trajectory[time][i_type][i_layer]
                    if len(startlayer) > 0:
                        tmp_config_hexatic = []
                        n_ordered = 0
                        # loop over particles
                        for i_particle in range(len(startlayer)):
                            if not self.particles[i_type][i_layer][i_particle]:
                                continue
                            if i_type == 0:
                                hexatic = self.calculateHexatic(i_type, i_particle, i_layer, time)
                            else:
                                hexatic = self.calulateDelaunayShapes(i_type, i_particle, i_layer, time)
                            single_hexatic.append(hexatic)
                            tmp_config_hexatic.append(hexatic)
                            byparticle_hexatic[i_particle].append(hexatic)
                            if hexatic > 0.5:
                                n_ordered = n_ordered + 1
                        config_hexatic.append( np.mean( np.array(tmp_config_hexatic) ) )
                        fraction_ordered.append( float(n_ordered) / float(len(startlayer)))

        # compute histogram
        print len(single_hexatic)
        vals, bins = np.histogram( np.array(single_hexatic), 20, (0.0, 1.0))
        print n_ordered, float(n_ordered) / len(single_hexatic), np.mean(single_hexatic), bins[np.argsort(vals)[-2]], len(single_hexatic)
        
        # write to file
        ofilename = self.params.output_path+"/hexatic_type"+str(i_type)+"_start"+str(times[0])+"_end"+str(max(times))+"_stride"+str(stride)+".dat"
        print "writing to", ofilename
        ofilewriter = csv.writer(open(ofilename, "wb"), delimiter=' ')
        for i in range(len(vals)):
            ofilewriter.writerow([bins[i], vals[i] / float(sum(vals))])


    def PlotGofR(self, itype, do_time_average=True, cutoff_dist=60.0, particle_density=0.0, binwidth=1.0):
        """public: compute and print simple g(r) for center-to-center distances for particles of same type

        output file name: self.params.output_path+"/gofr1d_type"+str(type)+"_start"+str(times[0])+"_end"+str(maxtime)+"_stride"+str(stride)+".dat"

        output file format: \n
        col1 distance \n
        col2 g(r) at that distance \n
        col3 stddev at that distance \n
        @param self The object pointer
        @param itype Int for particle type
        @param do_time_average Bool for whether to compute for current time only (False) or average over current and future times (true), default True
        """
        # set up times and stride
        mintime = self.params.time
        times = [mintime]
        if do_time_average:
            times = self.params.getTimes()
        stride = self.params.stride
        
        print "plotting g(r) for type", itype, "with stride", stride, "and times", times

        # set up lists for single-time and time-average distances
        Lx = self.params.getWidth(mintime)
        Ly = self.params.getHeight(mintime)
        maxr = min(min(Lx, Ly) * 0.5, cutoff_dist)
        bins = [i*binwidth for i in range(int(maxr/binwidth))]
        g_counts = np.zeros((len(bins), len(times)))
        g_norm_counts = np.zeros((len(bins), len(times)))

        # loop over all start times
        nparticles = 0.0
        for itime, time in enumerate(times):
            if time % stride == 0:
                print "time is", time

                Lx = self.params.getWidth(time)
                Ly = self.params.getHeight(time)

                # do top layer only (if it exists)
                ilayer=-1
                if len(self.trajectory[time][itype][1]) > 0:
                    ilayer = 1
                elif len(self.trajectory[time][itype][0]) > 0:
                    ilayer = 0
                else:
                    print "skipping time", time, "for type", itype, "because it's empty:", [ [len(plist) for plist in llist] for llist in self.trajectory[time] ]
                    return
                startlayer = self.trajectory[time][itype][ilayer]

                n_particles = 0.0
                # loop over particle pairs
                for i_particle in range(len(startlayer)):
                    i_pos = startlayer[i_particle]
                    if i_pos:
                        # skip particle if tagged
                        if not self.particles[itype][ilayer][i_particle]._tagged:
                            n_particles += 1.0
                            # loop over *all* particles, not just ones with larger index or tagged ones
                            for j_particle in range(len(startlayer)):
                                j_pos = startlayer[j_particle]
                                if j_pos:
                                    # calculate PBC distance and add to histograms
                                    xdist = self.PBC_diff(i_pos.x, j_pos.x, Lx)
                                    ydist = self.PBC_diff(i_pos.y, j_pos.y, Ly)
                                    dist = sqrt(xdist*xdist + ydist*ydist)
                                    if dist < maxr:
                                        for ibin in reversed(range(len(bins))):
                                            if dist > bins[ibin]:
                                                g_counts[ibin][itime] += 1
                                                break

        # bin and normalize g(r) histograms
        area_in_bin = [pi*(bin+binwidth)*(bin+binwidth) - pi*bin*bin for bin in bins]
        if (particle_density == 0):
            bg_particles_per_area = n_particles / Lx / Ly
        else:
            bg_particles_per_area = particle_density
        for itime in range(len(times)):
            g_norm_counts[:,itime] =  g_counts[:,itime] / area_in_bin[:] / bg_particles_per_area / n_particles

        # write to file
        ofilename = self.params.output_path+"/gofr1d_type"+str(itype)+"_start"+str(times[0])+"_end"+str(times[-1])+"_stride"+str(stride)+".dat"
        print "writing to", ofilename
        ofilewriter = csv.writer(open(ofilename, "wb"), delimiter=' ')
        for ibin in range(len(bins)):
            ofilewriter.writerow([bins[ibin]+0.5*binwidth, np.mean(g_norm_counts[ibin,:]), np.std(g_norm_counts[ibin,:]), np.sum(g_counts[ibin,:]) ])

    def PlotGofTheta(self, type=1, i_layer=1, do_time_average=True):
        """public: compute and print spatially-resolved g(r) for center-to-center distances for particles of same type

        output file name: self.params.output_path+"/psspatial_type"+str(type)+"_start"+str(times[0])+"_end"+str(maxtime)+"_stride"+str(stride)+".dat"

        output file format: \n
        col1 distance perpendicular to particle axis
        col2 distance parallel to particle axis
        col3 g(r) at that position

        @param self The object pointer
        @param type Int for particle type, default 1 (PSII)
        @param i_layer Int for layer id, default 1 (top, lumen-side-up)
        @param do_time_average Bool for whether to compute for current time only (False) or average over current and future times (true), default True
        """
        # spatially resolved orientation correlation
        # cf Bates & Frenkel 2000
        # g_l(r) = < cos ( l*(theta(0) - theta(r) ) >
        # for r only perpendicular or parallel to theta(0)
        # rank = l = 2, 4

        # set up times and stride
        mintime = self.params.time
        times = [mintime]
        if do_time_average:
            times = self.params.getTimes()
        maxtime = np.max(times)
        stride = self.params.stride
        
        # set up lists for single-time and time-average distances
        # only use with constant box size!
        Lx = self.params.getWidth(mintime)
        Ly = self.params.getHeight(mintime)
        maxr = min(firstLx, firstLy)  / 3.0
        bins = range(0, int(maxr))
        emptylayer = []
        for bin in bins: # rperp bins
            binlist = []
            for bin in bins: # rpar bins
                binlist.append([])
            emptylayer.append(binlist)
        timeav_g2 = emptylayer
        single_g2 = []
        timeav_gr = []
        for bin in bins: # rperp bins
            binlist = []
            for bin in bins: # rpar bins
                binlist.append(0)
            timeav_gr.append(binlist)

        # loop over all start times
        nparticles = 0.0
        for time in times:
            if time % stride == 0:
                print "time is", time
                single_g2 = emptylayer
                single_gr = emptylayer
                startlayer = self.trajectory[time][type][i_layer]
                if len(startlayer) > 0:
                    nparticles = 0.0
                    # loop over reference particles
                    for i_particle in range(len(startlayer)):
                        i_pos = startlayer[i_particle]
                        if i_pos:
                            # set up unit vectors parallel and perpendicular to reference particle
                            parvector = [ cos(i_pos.theta), sin(i_pos.theta) ]
                            perpvector = [ -sin(i_pos.theta), cos(i_pos.theta) ]
                            nparticles = nparticles + 1.0
                            # loop over other particles
                            for j_particle in range(i_particle+1, len(startlayer)):
                                j_pos = startlayer[j_particle]
                                if j_pos:
                                    # project separation vector onto par. and perp. vectors
                                    xdist = self.PBC_diff(i_pos.x, j_pos.x, Lx)
                                    ydist = self.PBC_diff(i_pos.y, j_pos.y, Ly)
                                    rpar = abs(xdist*parvector[0] + ydist*parvector[1])
                                    rperp = abs(xdist*perpvector[0] + ydist*perpvector[1])
                                    if rperp < 5 and rpar < 15:
                                        print "particles ipos", i_pos, "and jpos", j_pos, "at rpar =", rpar, "and rperp =", rperp
                                            # calculate orientation correlation and add to histogram
                                        g2theta = cos( 2.0 * (i_pos.theta - j_pos.theta) )
                                        g4theta = cos( 4.0 * (i_pos.theta - j_pos.theta) )
                                        if rpar < int(maxr):
                                            if rperp < int(maxr):
                                                single_g2[int(rperp)][int(rpar)].append(g2theta)
                                                timeav_g2[int(rperp)][int(rpar)].append(g2theta)
                                                timeav_gr[int(rperp)][int(rpar)] = timeav_gr[int(rperp)][int(rpar)] + 1
        # average data in each bin
        single_g2_z = np.zeros( ( len(bins), len(bins) ) )
        timeav_g2_z = np.zeros( ( len(bins), len(bins) ) )
        timeav_gr_z = np.zeros( ( len(bins), len(bins) ) )
        for perpbin in range(len(single_g2)):
            for parbin in range(len(single_g2[perpbin])):
                single_g2_z[parbin, perpbin] = np.mean(np.array(single_g2[perpbin][parbin]))
                timeav_g2_z[parbin, perpbin] = np.mean(np.array(timeav_g2[perpbin][parbin]))
                timeav_gr_z[parbin, perpbin] = 2.0 * float(timeav_gr[perpbin][parbin]) / float(nparticles)

        # plot as in http://matplotlib.sourceforge.net/examples/mplot3d/surface3d_demo.html
        print "about to try to plot"
#        fig = plt.figure()
#        plt.suptitle("orientation correlation function starting at time "+str(times[0])+" for file "+self._parent.fileControlPanel.dir_text.GetLabel())
#        plt.subplot(121)
#        cax = plt.imshow(timeav_g2_z, cmap=cm.RdBu, origin='lower', vmin=-1, vmax=1)
#        plt.xlabel("r_perp in nm (rod width 12 nm)")
#        plt.ylabel("r_parallel in nm (rod length 26.5 nm)")
#        plt.title("g_2")
#        cbar = fig.colorbar(cax)
#        plt.grid()
#        plt.subplot(122)
#        cax = plt.imshow(timeav_gr_z, cmap=cm.RdBu, origin='lower', vmin=0, vmax=2)
#        plt.xlabel("r_perp in nm (rod width 12 nm)")
#        plt.ylabel("r_parallel in nm (rod length 26.5 nm)")
#        cbar = fig.colorbar(cax, extend='max')
#        plt.title("g(r)")
#        plt.grid()
#        plt.show()

        # write to file
        ofilename = self.params.output_path+"/psspatial_type"+str(type)+"_start"+str(times[0])+"_end"+str(maxtime)+"_stride"+str(stride)+".dat"
        print "writing to", ofilename
        ofilewriter = csv.writer(open(ofilename, "wb"), delimiter=' ')
        lastbin = len(bins) - 1
        for i in range(len(bins)):
            for j in range(len(bins)):
                ofilewriter.writerow([bins[i], bins[j], timeav_gr_z[i][j]
#                                  timeav_g2_z[0][i],
#                                  timeav_g2_z[i][0],
#                                  (timeav_g2_z[i][i] + timeav_g2_z[min([i+1,lastbin])][max([i-1,0])] +
#                                   timeav_g2_z[max([i-1,0])][min([i+1,lastbin])] +
#                                   timeav_g2_z[min([i+2,lastbin])][max([i-2,0])] +
#                                   timeav_g2_z[max([i-2,0])][min([i+2,lastbin])])/5.0,
#                                  timeav_gr_z[0][i],
#                                  timeav_gr_z[i][0],
#                                  (timeav_gr_z[i][i] + timeav_gr_z[min([i+1,lastbin])][max([i-1,0])] +
#                                   timeav_gr_z[max([i-1,0])][min([i+1,lastbin])] +
#                                   timeav_gr_z[min([i+2,lastbin])][max([i-2,0])] +
#                                   timeav_gr_z[max([i-2,0])][min([i+2,lastbin])])/5.0
                                     ])

    def PlotMsd(self, type):
        """public: compute and print mean square displacement for all particles of given type
        
        output file name: self.params.output_path+"/rsq_type"+str(type)+"_start"+str(times[0])+"_end"+str(maxtime)+"_stride"+str(stride)+".dat"

        output file format: \n
        col1 time \n
        col2 MSD for all particles \n
        col3 MSD for phosphorylated particles in grana \n
        col4 MSD for unphosphorylated particles in grana \n
        col5 MSD for all particles in stroma \n
        col6 non-Gaussian parameter alpha of squared-displacement distribution, 1=Gaussian \n

        deprecated output file: histogram of displacements at shortest time

        @param self The object pointer
        @param type Int for particle type
        """
        # set up times and stride
        mintime = self.params.time
        times = self.params.getTimes()
        maxtime = np.max(times)
        stride = self.params.stride

        # shouldn't be computing diffusion with changing box size
        # so ok to get box size once
        Lx = self.params.getWidth(mintime)
        Ly = self.params.getHeight(mintime)

        print "plotting msd for type", type, "with stride", stride, "and times", times

        # set up r^2 and r^4 lists for MSD
        sum_r_sq = [ 0.0 ]
        count_r_sq = [ 0.0 ]
        sum_r_4th = [ 0.0 ]
        count_r_4th = [ 0.0 ]
        sum_grana_unphos_r_sq = [ 0.0 ]
        count_grana_unphos_r_sq = [ 0.0 ]
        sum_grana_phos_r_sq = [ 0.0 ]
        count_grana_phos_r_sq = [ 0.0 ]
        sum_stroma_r_sq = [ 0.0 ]
        count_stroma_r_sq = [ 0.0 ]

        # set up r^2 and r^4 lists for MME (Tejedor , ..., Metzler 2010)
#        max_r_sq = [ [], [] ]
#        max_r_4th = [ [], [] ]

        # loop over all start times
        for starttime in times:
            if starttime % stride == 0:
                print "start time is", starttime #, len(sum_r_sq), sum_r_sq[0], sum_r_sq[1], count_r_sq[1]
                # loop over layers
                for i_layer in range(len(self.trajectory[starttime][type])):
                    startlayer = self.trajectory[starttime][type][i_layer]
                    # loop over particles
                    for i_particle in range(len(startlayer)):
                        # get position from trajectory
                        startpos = startlayer[i_particle]
                        pbc_image = [0, 0]
                        prev_corrpos = startpos
                        max_rsq = 0
                        particle = self.particles[type][i_layer][i_particle]
                        if startpos:
                            # loop over later times (corrtime)
                            i_corr = 1
                            while starttime + i_corr*stride <= maxtime:
                                # set up later time
                                corrtime = starttime + i_corr*stride
                                if i_corr >= len(sum_r_sq):
                                    sum_r_sq.append(0.0)
                                    count_r_sq.append(0.0)
                                    sum_grana_unphos_r_sq.append(0.0)
                                    count_grana_unphos_r_sq.append(0.0)
                                    sum_grana_phos_r_sq.append(0.0)
                                    count_grana_phos_r_sq.append(0.0)
                                    sum_stroma_r_sq.append(0.0)
                                    count_stroma_r_sq.append(0.0)
                                    sum_r_4th.append(0.0)
                                    count_r_4th.append(0.0)

                                # get position after lag
                                if corrtime in self.trajectory:
                                    corrpos = self.trajectory[corrtime][type][i_layer][i_particle]
                                    pbc_image = corrpos.update_images(pbc_image, prev_corrpos, Lx, Ly)

                                    # get rsq using appropriate PBC positions
                                    #                                xdist = self.PBC_diff(startpos[0], corrpos[0], Lx) + pbc_image[0] * Lx
#                                ydist = self.PBC_diff(startpos[1], corrpos[1], Ly) + pbc_image[1] * Ly
                                    xdist = startpos.x - (corrpos.x + pbc_image[0] * Lx)
                                    ydist = startpos.y - (corrpos.y + pbc_image[1] * Ly)
                                    rsq = (xdist*xdist + ydist*ydist) * 0.000001
                                
                                    # test if rsq is maximal
                                    if rsq > max_rsq:
                                        max_rsq = rsq

                                    # add rsq and max_rsq to lists
                                    sum_r_sq[i_corr] += rsq
                                    count_r_sq[i_corr] += 1
                                    if startpos.region and corrpos.region:
                                        if startpos.region==1 and corrpos.region==1:
                                            if particle._phos:
                                                sum_grana_phos_r_sq[i_corr] += rsq
                                                count_grana_phos_r_sq[i_corr] += 1
                                            else:
                                                sum_grana_unphos_r_sq[i_corr] += rsq
                                                count_grana_unphos_r_sq[i_corr] += 1
                                        elif startpos.region==2 and corrpos.region==2:
                                            sum_stroma_r_sq[i_corr] += rsq
                                            count_stroma_r_sq[i_corr] += 1
                                    sum_r_4th[i_corr] += rsq*rsq
                                    count_r_4th[i_corr] += 1
                                else:
                                    print "missing corrtime", corrtime, " so skipping it"

                                # prepare for next loop
                                i_corr = i_corr+1
                                prev_corrpos = corrpos

                                    

        # make mean and percentile lists
        xaxis_times = [ i_t * stride * 0.000001 for i_t in range(len(sum_r_sq)) ]
        means_rsq = [ float(sum_r_sq[i_t]) / max(1.0, float(count_r_sq[i_t])) for i_t in range(len(sum_r_sq)) ]
        grana_phos_means_rsq = [ float(sum_grana_phos_r_sq[i_t]) /  max(1.0, float(count_grana_phos_r_sq[i_t])) for i_t in range(len(sum_r_sq)) ]
        grana_unphos_means_rsq = [ float(sum_grana_unphos_r_sq[i_t]) /  max(1.0, float(count_grana_unphos_r_sq[i_t])) for i_t in range(len(sum_r_sq)) ]
        stroma_means_rsq = [ float(sum_stroma_r_sq[i_t]) /  max(1.0, float(count_stroma_r_sq[i_t])) for i_t in range(len(sum_r_sq)) ]
        alphas_r4th = [ float(sum_r_4th[i_t]) /  max(1e-7, float(count_r_4th[i_t])) / max(1e-7, means_rsq[i_t]* means_rsq[i_t]) / 3.0 - 1.0   for i_t in range(len(sum_r_sq)) ]

        # write to file
        ofilename = self.params.output_path+"/rsq_type"+str(type)+"_start"+str(times[0])+"_end"+str(maxtime)+"_stride"+str(stride)+".dat"
        print "writing to", ofilename
        ofilewriter = csv.writer(open(ofilename, "wb"), delimiter=' ')
        for i in range(1, len(xaxis_times)):
            ofilewriter.writerow([xaxis_times[i], means_rsq[i] , grana_phos_means_rsq[i], grana_unphos_means_rsq[i], stroma_means_rsq[i], alphas_r4th[i]])

      #  ofilename2 = self._parent.fileControlPanel.output_path+"/disp_hist_type"+str(type)+"_start"+str(times[0])+"_end"+str(maxtime)+"_stride"+str(stride)+".dat"
      #  print "writing to", ofilename2
      #  ofilewriter2 = csv.writer(open(ofilename2, "wb"), delimiter=' ')
      #  for i in range(len(bins)-1):
      #      ofilewriter2.writerow([ bins[i], hist_all[i], hist_grana_phos[i], hist_grana_unphos[i], hist_stroma[i] ])


    def PlotSingleMsd(self, particle):
        """public: compute and return Numpy array with mean square displacement for a single particle
        
        @param self The object pointer
        @param particle Particle object
        @retval nparray Numpy array [times, msds, non-Gaussian parameters]
        """
        # set up times and stride
        mintime = self.params.time
        times = self.params.getTimes()
        maxtime = np.max(times)
        stride = self.params.stride

        # shouldn't be computing diffusion with changing box size
        # so ok to get box size once
        Lx = self.params.getWidth(mintime)
        Ly = self.params.getHeight(mintime)

        print "plotting msd for type", particle._type, "with stride", stride, "and times", times

        # set up r^2 and r^4 lists for MSD
        sum_r_sq = [ 0.0 ]
        count_r_sq = [ 0.0 ]
        sum_r_4th = [ 0.0 ]
        count_r_4th = [ 0.0 ]

        # loop over all start times
        for starttime in times:
            if starttime % stride == 0:
                print "start time is", starttime #, len(sum_r_sq), sum_r_sq[0], sum_r_sq[1], count_r_sq[1]
                # get position from trajectory
                startpos = self.trajectory[starttime][particle._type][particle._layer][particle._id]
                pbc_image = [0, 0]
                prev_corrpos = startpos
                max_rsq = 0
                if startpos:
                    # loop over later times (corrtime)
                    i_corr = 1
                    while starttime + i_corr*stride <= maxtime:
                        # set up later time
                        corrtime = starttime + i_corr*stride
                        if i_corr >= len(sum_r_sq):
                            sum_r_sq.append(0.0)
                            count_r_sq.append(0.0)
                            sum_r_4th.append(0.0)
                            count_r_4th.append(0.0)

                        # get position after lag
                        if corrtime in self.trajectory:
                            corrpos = self.trajectory[corrtime][particle._type][particle._layer][particle._id]
                            pbc_image = corrpos.update_images(pbc_image, prev_corrpos, Lx, Ly)
                            
                            # get rsq using appropriate PBC positions
                            #                                xdist = self.PBC_diff(startpos[0], corrpos[0], Lx) + pbc_image[0] * Lx
                            #                                ydist = self.PBC_diff(startpos[1], corrpos[1], Ly) + pbc_image[1] * Ly
                            xdist = startpos.x - (corrpos.x + pbc_image[0] * Lx)
                            ydist = startpos.y - (corrpos.y + pbc_image[1] * Ly)
                            rsq = (xdist*xdist + ydist*ydist) * 0.000001
                            
                            # test if rsq is maximal
                            if rsq > max_rsq:
                                max_rsq = rsq

                            # add rsq and max_rsq to lists
                            sum_r_sq[i_corr] += rsq
                            count_r_sq[i_corr] += 1
                            sum_r_4th[i_corr] += rsq*rsq
                            count_r_4th[i_corr] += 1
                        else:
                            print "missing corrtime", corrtime, " so skipping it"
    
                        # prepare for next loop
                        i_corr = i_corr+1
                        prev_corrpos = corrpos

        # make mean and percentile lists
        xaxis_times = [ i_t * stride * 0.000001 for i_t in range(len(sum_r_sq)) ]
        means_rsq = [ float(sum_r_sq[i_t]) / max(1.0, float(count_r_sq[i_t])) for i_t in range(len(sum_r_sq)) ]
        alphas_r4th = [ float(sum_r_4th[i_t]) /  max(1e-7, float(count_r_4th[i_t])) / max(1e-7, means_rsq[i_t]* means_rsq[i_t]) / 3.0 - 1.0   for i_t in range(len(sum_r_sq)) ]

        return np.array([xaxis_times, means_rsq, alphas_r4th])

    def PlotAngleCorr(self, itype=1, dbinwidth=1.0):
        """public: compute histogram of orientation correlation (ie histogram of cos(theta)) over distance

        uncorrelated value (long-distance limit) is 2/pi

        output file name: self.params.output_path+"/anglecorr_type"+str(itype)+"_start"+str(times[0])+"_end"+str(maxtime)+"_stride"+str(stride)+".dat"

        output file format: \n
        col1 bin in nm \n
        col2 angle-angle correlation for in-plane \n
        col3 stddev \n
        col4 counts \n

        @param self The object pointer
        @param itype Int for particle type
        @param dbinwidth Float for width of distance bin, in nm
        """
        # set up times and stride
        mintime = self.params.time
        times = self.params.getTimes()
        maxtime = np.max(times)
        stride = self.params.stride

        print "plotting angle-angle correlation over distance for type", itype, "with stride", stride, "and times", times

        # set up arrays for bins and counts
        bins = [i*dbinwidth for i in range(int(100/dbinwidth))]
        corrs = np.zeros((len(bins), len(times)))
        counts = np.zeros((len(bins), len(times)))
        norms = np.zeros((len(bins), len(times)))

        # loop over particles in one layer
        for itime, starttime in enumerate(times):
            if starttime % stride == 0:
                print "time is", starttime

                # choose reference layer (top if possible)
                if len(self.trajectory[starttime][itype][1]) > 0:
                    ilayer = 1
                else:
                    ilayer = 0

                # set up positions list
                positions = []
                for i_particle in range(len(self.trajectory[starttime][itype][ilayer])):
                    pos = self.trajectory[starttime][itype][ilayer][i_particle]
                    if pos:
                        positions.append( (itype, i_particle, self.particles[itype][ilayer][i_particle]._tagged, pos) )

                # loop over particles
                for (ref_itype, ref_ip, ref_tagged, ref_pos) in positions:

                    # skip particle if tagged
                    if not ref_tagged:
                            
                        # loop over other particles
                        for (other_itype, other_ip, other_tagged, other_pos) in positions:
                            if other_itype == ref_itype and other_ip == ref_ip:
                                pass
                            else:
                                xdist = self.PBC_diff(ref_pos.x, other_pos.x, self.params.getWidth(starttime))
                                ydist = self.PBC_diff(ref_pos.y, other_pos.y, self.params.getHeight(starttime))
                                r = sqrt(xdist*xdist + ydist*ydist)
                                costheta = abs(cos(ref_pos.theta - other_pos.theta))                                    

                                # add cos(theta) to histograms
                                for ibin in reversed(range(len(bins))):
                                    if r > bins[ibin]:
                                        counts[ibin][itime] += 1
                                        corrs[ibin][itime] += costheta
                                        break

        # normalize correlation
        for itime in range(len(times)):
            for ibin in range(len(bins)):
                norms[ibin,itime] = corrs[ibin,itime] / max(counts[ibin,itime], 1.0)

        # print to files
        distfilename = self.params.output_path+"/anglecorr_type"+str(itype)+"_start"+str(times[0])+"_end"+str(maxtime)+"_stride"+str(stride)+".dat"
        distfilewriter = csv.writer(open(distfilename, "wb"), delimiter=' ')
        for ibin in range(len(bins)):
            distfilewriter.writerow([bins[ibin]+0.5*dbinwidth,
                                     np.mean(norms[ibin]), np.std(norms[ibin]), np.sum(counts[ibin]),
                                     ])
                                     
    def PlotNeighborsWithin(self, itype, rcut, dbinwidth=1.0):
        """public: compute histogram of onumber of particles within rcut nm of a central particle

        uncorrelated value (long-distance limit) is 2/pi

        output file name: self.params.output_path+"/nnwithin_type"+str(itype)+"_start"+str(times[0])+"_end"+str(maxtime)+"_stride"+str(stride)+".dat"

        output file format: \n
        col1 bin in number of particles \n
        col2 prob for number of neighbors \n
        col3 stddev \n
        col4 counts \n

        @param self The object pointer
        @param itype Int for particle type
        @param dbinwidth Float for width of distance bin, in nm
        """
        # set up times and stride
        mintime = self.params.time
        times = self.params.getTimes()
        maxtime = np.max(times)
        stride = self.params.stride

        print "plotting number of neighbors within", rcut, "for type", itype, "with stride", stride, "and times", times

        # set up arrays for bins and counts
        bins = [i*dbinwidth for i in range(int(35/dbinwidth))]
        counts = np.zeros((len(bins), len(times)))
        norms = np.zeros((len(bins), len(times)))

        # loop over particles in one layer
        for itime, starttime in enumerate(times):
            if starttime % stride == 0:
                print "time is", starttime

                # choose reference layer (top if possible)
                if len(self.trajectory[starttime][itype][1]) > 0:
                    ilayer = 1
                else:
                    ilayer = 0

                # set up positions list
                positions = []
                for i_particle in range(len(self.trajectory[starttime][itype][ilayer])):
                    pos = self.trajectory[starttime][itype][ilayer][i_particle]
                    if pos:
                        positions.append( (itype, i_particle, self.particles[itype][ilayer][i_particle]._tagged, pos) )
                # add bound LHC if target is free LHC
                if itype == 0:
                    for i_particle in range(len(self.trajectory[starttime][2][ilayer])):
                        pos = self.trajectory[starttime][2][ilayer][i_particle]
                        if pos:
                            positions.append( (2, i_particle, self.particles[2][ilayer][i_particle]._tagged, pos) )

                # loop over particles
                for (ref_itype, ref_ip, ref_tagged, ref_pos) in positions:

                    # skip particle if tagged
                    if not ref_tagged:
                        # set up counter for number of neighbors
                        nn = 0
                        
                        # loop over other particles
                        for (other_itype, other_ip, other_tagged, other_pos) in positions:
                            if other_itype == ref_itype and other_ip == ref_ip:
                                pass
                            else:
                                xdist = self.PBC_diff(ref_pos.x, other_pos.x, self.params.getWidth(starttime))
                                ydist = self.PBC_diff(ref_pos.y, other_pos.y, self.params.getHeight(starttime))
                                r = sqrt(xdist*xdist + ydist*ydist)
                                # increment counter
                                if r < rcut:                          
                                    nn += 1
                                    
                        # add counter to histograms
                        for ibin in reversed(range(len(bins))):
                            if nn > bins[ibin]:
                                counts[ibin][itime] += 1
                                break

        # normalize correlation
        for itime in range(len(times)):
            for ibin in range(len(bins)):
                norms[ibin,itime] = counts[ibin,itime] / max(sum(counts[:,itime]), 1.0)

        # print to files
        distfilename = self.params.output_path+"/nnwithin_type"+str(itype)+"_start"+str(times[0])+"_end"+str(maxtime)+"_stride"+str(stride)+".dat"
        distfilewriter = csv.writer(open(distfilename, "wb"), delimiter=' ')
        for ibin in range(len(bins)):
            distfilewriter.writerow([bins[ibin]+0.5*dbinwidth,
                                     np.mean(norms[ibin]), np.std(norms[ibin]), np.sum(counts[ibin]),
                                     ])



    def PlotNearestNbhrs(self, ptype, dbinwidth=1.0, tbinwidth=5.0):
        """public: compute nearest-neighbor distance and angle correlations

        theta output file name: self.params.output_path+"/nntheta_type"+str(ptype)+"_start"+str(times[0])+"_end"+str(maxtime)+"_stride"+str(stride)+".dat"

        theta file format: \n
        col1 bin in radians \n
        col2 probability of NN angle difference for in-plane neighbors \n
        col3 stddev for in-plane \n
        col4 count for in-plane \n
        col5 probability of NN angle difference for stacked neighbors \n
        col6 stddev for stacked \n
        col7 count for stacked \n

        dist output file name: self.params.output_path+"/nndist_type+"str(ptype)+"_start"+str(times[0])+"_end"+str(maxtime)+"_stride"+str(stride)+".dat"

        dist file format: \n
        col1 bin in nm \n
        col2 probability of NN distance for in-plane neighbors \n
        col3 stddev for in-plane \n
        col4 count for in-plane \n
        col5 probability of NN distance for stacked neighbors \n
        col6 stddev for stacked \n
        col7 count for stacked \n

        @param self The object pointer
        @param ptype Int for particle type
        @param dbinwidth Float for width of distance bin, in nm
        @param tbinwidth Float for width of angle bin, in degrees
        """
        # set up times and stride
        mintime = self.params.time
        times = self.params.getTimes()
        maxtime = np.max(times)
        stride = self.params.stride

        print "plotting NN distance and orientations for type", ptype, "with stride", stride, "and times", times

        # set up arrays for bins and counts
        degree_bins = [i*tbinwidth for i in range(int(90.0/tbinwidth))]
        theta_bins = [cos(deg/180.0*pi) for deg in degree_bins]
        stacked_theta_counts = np.zeros((len(theta_bins), len(times)))
        inplane_theta_counts = np.zeros((len(theta_bins), len(times)))
        stacked_theta_norms = np.zeros((len(theta_bins), len(times)))
        inplane_theta_norms = np.zeros((len(theta_bins), len(times)))
        dist_bins = [i*dbinwidth for i in range(int(50/dbinwidth))]
        stacked_dist_counts = np.zeros((len(dist_bins), len(times)))
        inplane_dist_counts = np.zeros((len(dist_bins), len(times)))
        stacked_dist_norms = np.zeros((len(dist_bins), len(times)))
        inplane_dist_norms = np.zeros((len(dist_bins), len(times)))

        # loop over particles in one layer
        for itime, starttime in enumerate(times):
            if starttime % stride == 0:
                print "time is", starttime

                # choose reference layer (top if possible)
                if len(self.trajectory[starttime][ptype][1]) > 0:
                    istartlayer = 1
                    iotherlayer = 0
                else:
                    istartlayer = 0
                    iotherlayer = 1

                # correlate both free and embedded LHCII
                types = [ptype]
                if ptype == 0:
                    types.append(2)

                # set up position lists
                positions = [ [], [] ]
                for itype in types:
                    for ilayer in [istartlayer, iotherlayer]:
                        for i_particle in range(len(self.trajectory[starttime][itype][ilayer])):
                            pos = self.trajectory[starttime][itype][ilayer][i_particle]
                            if pos:
                                positions[ilayer].append( (itype, i_particle, self.particles[itype][ilayer][i_particle]._tagged, pos) )

                # loop over particles
                for (ref_itype, ref_ip, ref_tagged, ref_pos) in positions[istartlayer]:

                    # skip particle if tagged
                    if not ref_tagged:
                            
                        # find in-plane neighbor
                        minrsq = 70*70
                        theta = -100
                        for (other_itype, other_ip, other_tagged, other_pos) in positions[istartlayer]:
                            if other_itype == ref_itype and other_ip == ref_ip:
                                pass
                            else:
                                xdist = self.PBC_diff(ref_pos.x, other_pos.x, self.params.getWidth(starttime))
                                ydist = self.PBC_diff(ref_pos.y, other_pos.y, self.params.getHeight(starttime))
                                rsq = xdist*xdist + ydist*ydist
                                if rsq < minrsq:
                                    minrsq = rsq
                                    theta = abs(cos(ref_pos.theta - other_pos.theta))                                    

                        # add to histograms
                        for ibin in reversed(range(len(theta_bins))):
                            if theta < theta_bins[ibin]:
                                inplane_theta_counts[ibin][itime] += 1
                                break
                        for ibin in reversed(range(len(dist_bins))):
                            if sqrt(minrsq) > dist_bins[ibin]:
                                inplane_dist_counts[ibin][itime] += 1
                                break

                        # find stacking neighbor
                        minrsq = 70*70
                        theta = -100
                        for (other_itype, other_ip, other_tagged, other_pos) in positions[iotherlayer]:
                            xdist = self.PBC_diff(ref_pos.x, other_pos.x, self.params.getWidth(starttime))
                            ydist = self.PBC_diff(ref_pos.y, other_pos.y, self.params.getHeight(starttime))
                            rsq = xdist*xdist + ydist*ydist
                            if rsq < minrsq:
                                minrsq = rsq
                                theta = abs(cos(ref_pos.theta - other_pos.theta))                                    
                                        
                        # add to histograms
                        for ibin in reversed(range(len(theta_bins))):
                            if theta < theta_bins[ibin]:
                                stacked_theta_counts[ibin][itime] += 1
                                break
                        for ibin in reversed(range(len(dist_bins))):
                            if sqrt(minrsq) > dist_bins[ibin]:
                                stacked_dist_counts[ibin][itime] += 1
                                break
                         
        # normalize counts
        for itime in range(len(times)):
            inplane_theta_norms[:,itime] = inplane_theta_counts[:,itime] / max(np.sum(inplane_theta_counts[:,itime]), 1.0)
            stacked_theta_norms[:,itime] = stacked_theta_counts[:,itime] / max(np.sum(stacked_theta_counts[:,itime]), 1.0)
            inplane_dist_norms[:,itime] = inplane_dist_counts[:,itime] / max(np.sum(inplane_dist_counts[:,itime]), 1.0)
            stacked_dist_norms[:,itime] = stacked_dist_counts[:,itime] / max(np.sum(stacked_dist_counts[:,itime]), 1.0)

        # print to files
        if ptype==1:
            thetafilename = self.params.output_path+"/nntheta_type"+str(ptype)+"_start"+str(times[0])+"_end"+str(maxtime)+"_stride"+str(stride)+".dat"
            thetafilewriter = csv.writer(open(thetafilename, "wb"), delimiter=' ')
            for ibin in range(len(theta_bins)):
                thetafilewriter.writerow([degree_bins[ibin]+0.5*tbinwidth,
                                          np.mean(inplane_theta_norms[ibin]), np.std(inplane_theta_norms[ibin]), np.sum(inplane_theta_counts[ibin]),
                                          np.mean(stacked_theta_norms[ibin]), np.std(stacked_theta_norms[ibin]), np.sum(stacked_theta_counts[ibin])
                                          ])
                                 
        distfilename = self.params.output_path+"/nndist_type"+str(ptype)+"_start"+str(times[0])+"_end"+str(maxtime)+"_stride"+str(stride)+".dat"
        distfilewriter = csv.writer(open(distfilename, "wb"), delimiter=' ')
        for ibin in range(len(dist_bins)):
            distfilewriter.writerow([dist_bins[ibin]+0.5*dbinwidth,
                                     np.mean(inplane_dist_norms[ibin]), np.std(inplane_dist_norms[ibin]), np.sum(inplane_dist_counts[ibin]),
                                     np.mean(stacked_dist_norms[ibin]), np.std(stacked_dist_norms[ibin]), np.sum(stacked_dist_counts[ibin])
                                     ])

        # print summary to screen
        probs = np.mean(inplane_dist_norms, axis=1)
        print "most probable in-plane nndist and prob:", dist_bins[np.argmax(probs)]+0.5*dbinwidth, np.max(probs)
        print "average in-plane nndist:", np.average(dist_bins, weights=probs)+0.5*dbinwidth

    def PlotNematic(self, type):
        """deprecated: compute and print histogram of local nematic order parameters
        @param self The object pointer
        @param type Int for particle type
        """
        # set up times and stride
        mintime = self.params.time
        times = self.params.getTimes()
        maxtime = np.max(times)
        stride = self.params.stride

        print "plotting nematic for type", type, "with stride", stride, "and times", times

        # set up arrays for raw theta, theta diff to stacked neighbor, and theta diff to in-plane neighbor
        nematics = []
        nbhrs = []

        # loop over times
        for time in times:
            if time % stride == 0:
                print "time is", time
                # loop over layers
                for i_layer in range(len(self.trajectory[time][type])):
                    # loop over particles
                    for i_particle in range(len(self.trajectory[time][type][i_layer])):
                        pos = self.trajectory[time][type][i_layer][i_particle]
                        if pos:
                            nematic, nbhr = self.calculateNematic(type, i_particle, i_layer, time, True)
                            nematics.append(nematic)
                            nbhrs.append(nbhr)

        # calculate 2D histogram of nematic and # nbhrs
        hist, xedges, yedges = np.histogram2d(nematics, nbhrs, bins=[50, max(nbhrs)-min(nbhrs)])
        print hist.shape, xedges.shape, yedges.shape

        # output
        filename = self.params.output_path+"/nematic_hist_start"+str(times[0])+"_end"+str(maxtime)+"_stride"+str(stride)+".dat"
        filewriter = csv.writer(open(filename, "wb"), delimiter=' ')
        for ix in range(len(xedges)-1):
            for iy in range(len(yedges)-1):
                filewriter.writerow([xedges[ix], yedges[iy+1], hist[ix][iy]])
            filewriter.writerow([])     

    def PlotRegions(self):
        """public: computes and prints info about which particles are in which regions

        traj output file name: self.params.output_path+"/regions_traj_start"+str(times[0])+"_end"+str(maxtime)+"_stride"+str(stride)+".dat"

        traj file format: \n
        col1 time \n
        col2 number of unphosphorylated LHCII in grana \n
        col3 number of PSII in grana \n
        col4 number of phosphorylated LHCII in grana \n
        col5 number of unphosphorylated LHCII in stroma \n
        col6 number of PSII in stroma \n
        col7 number of phosphorylated LHCII in stroma \n        
        col8 number of unphosphorylated LHCII in disallowed region \n
        col9 number of PSII in disallowed region \n
        col10 number of phosphorylated LHCII in disallowed region \n       

        stats output file name: self.params.output_path+"/regions_stats_start"+str(times[0])+"_end"+str(maxtime)+"_stride"+str(stride)+".dat"

        stats file format: \n
        col1 type (loop over 0,1) \n
        col2 region (loop over 0,1,2) \n
        col3 average number of those particles in that region \n
        col4 standard deviation of number of those particles in that region \n
        @param self The object pointer
        """
        # set up times and stride
        mintime = self.params.time
        times = self.params.getTimes()
        maxtime = np.max(times)
        stride = self.params.stride

        print "plotting regions with stride", stride, "and times", times

        regions = [0, 1, 2]
        types = [0, 1]
        phostypes = [0, 1, 2] # type '2' here is phosphorylated free LHCII
        data = np.zeros( ( len(phostypes), len(regions), len(times) ) )

        # loop over times
        for i_time in range(len(times)):
            starttime = times[i_time]
            if starttime % stride == 0:
                print "time is", starttime
                # loop over types
                for type in types:
                    # loop over layers
                    for i_layer in range(len(self.trajectory[starttime][type])):
                        startlayer = self.trajectory[starttime][type][i_layer]
                        # loop over particles
                        for i_particle in range(len(startlayer)):
                            # get position from trajectory
                            startpos = startlayer[i_particle]
                            if startpos:
                                phostype = type
                                if self.particles[type][i_layer][i_particle]._phos:
                                    phostype = 2
                                if startpos.region:
                                    region = startpos.region
                                    data[phostype][region][i_time] += 1
                                else:
                                    print type, i_layer, startpos, "doesn't include region"

        
        # average for each region for each type
        means = data.mean(axis=2)
        stddevs = data.std(axis=2)

        # write to file
        trajfilename = self.params.output_path+"/regions_traj_start"+str(times[0])+"_end"+str(maxtime)+"_stride"+str(stride)+".dat"
        print "writing traj to", trajfilename
        trajfilewriter = csv.writer(open(trajfilename, "wb"), delimiter=' ')
        for i_time in range(len(times)):
            row = [times[i_time]]
            for type in phostypes:
                for region in regions:
                    row.append(data[type][region][i_time])
            trajfilewriter.writerow(row)

        statsfilename = self.params.output_path+"/regions_stats_start"+str(times[0])+"_end"+str(maxtime)+"_stride"+str(stride)+".dat"
        print "writing stats to", statsfilename
        statsfilewriter = csv.writer(open(statsfilename, "wb"), delimiter=' ')
        row = []
        for type in phostypes:
            for region in regions:
                statsfilewriter.writerow([ type, region, means[type][region], stddevs[type][region] ])

    def PlotClusters(self):
        """public: compute and print info about energy transfer clusters

        traj output file name: self.params.output_path+"/cluster_traj_start"+str(times[0])+"_end"+str(maxtime)+"_stride"+str(stride)+".dat"

        traj output file format: \n
        col1 time \n
        col2 layer id \n
        col3 number of LHCII in largest cluster \n
        col4 number of PSII in largest cluster \n
        
        antenna output file name: self.params.output_path+"/antenna_hist_start"+str(times[0])+"_end"+str(maxtime)+"_stride"+str(stride)+".dat"

        antenna output file format: \n
        col1 number of Chl per PSII reaction center \n
        col2 number of reaction centers with that many Chl \n
        col3 fraction of reaction centers with that many Chl \n

        cluster output file name: self.params.output_path+"/cluster_hist_start"+str(times[0])+"_end"+str(maxtime)+"_stride"+str(stride)+".dat"

        cluster output file format: \n
        col1 number of LHCII in cluster \n
        col2 number of PSII in cluster \n
        col3 number of clusters with that many particles

        @param self The object pointer
        """
        # set up times and stride
        mintime = self.params.time
        times = self.params.getTimes()
        maxtime = np.max(times)
        stride = self.params.stride

        print "plotting regions with stride", stride, "and times", times

        # set up Chl params and histograms
        types = [0, 1]
        chl_per_type = [42.0, 200.0]
        antenna_bin_min = 0.5*chl_per_type[1]
        antenna_bin_width = 0.2*chl_per_type[0]
        antenna_size_bins = [antenna_bin_min+i*antenna_bin_width for i in range(70)]
        antenna_size_hist = [0.0 for bin in antenna_size_bins]
        cluster_1d_bins = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 800, 850, 900, 950, 1000, 1200, 1400, 1600, 1800, 2000]
        cluster_2d_hist = [ [0 for bin in cluster_1d_bins] for bin in cluster_1d_bins ]

        # set up traj file
        trajfilename = self.params.output_path+"/cluster_traj_start"+str(times[0])+"_end"+str(maxtime)+"_stride"+str(stride)+".dat"
        print "writing largest cluster trajectory to", trajfilename
        trajfilewriter = csv.writer(open(trajfilename, "wb"), delimiter=' ')

        # loop over times
        for i_time in range(len(times)):
            starttime = times[i_time]
            if starttime % stride == 0:
                print "time is", starttime
                # loop over layers
                for i_layer in range(len(self.trajectory[starttime][0])):
                    # build clusters, if you haven't already
                    if starttime not in self.cluster_sizes:
                        self.buildETClusters(starttime, i_layer)
                    if i_layer not in self.cluster_sizes[starttime]:
                        self.buildETClusters(starttime, i_layer)
                    for cluster in self.cluster_sizes[starttime][i_layer]:
                        # build per-ps antenna size histogram
                        if cluster[1] > 0:
                            full_antenna_size = (chl_per_type[0]*cluster[0] + chl_per_type[1]*cluster[1])
                            antenna_size_per_ps = full_antenna_size / (2.0*cluster[1])
                            bin = int((antenna_size_per_ps-antenna_bin_min)/antenna_bin_width)
                            antenna_size_hist[bin] += 2.0*cluster[1]
                        
                        # build cluster size histogram
                        ps_bin = 0
                        lhc_bin = 0
                        while cluster_1d_bins[lhc_bin] < cluster[0]:
                            if lhc_bin < len(cluster_1d_bins):
                                lhc_bin += 1
                            else:
                                break
                        while cluster_1d_bins[ps_bin] < cluster[1]:
                            if ps_bin < len(cluster_1d_bins):
                                ps_bin += 1
                            else:
                                break
                        cluster_2d_hist[lhc_bin][ps_bin] += 1
                    # write largest cluster to file
                    trajfilewriter.writerow([starttime, i_layer]+self.cluster_sizes[starttime][i_layer][0])

        # write data
        antennafilename = self.params.output_path+"/antenna_hist_start"+str(times[0])+"_end"+str(maxtime)+"_stride"+str(stride)+".dat"
        print "writing antenna hist to", antennafilename
        antennafilewriter = csv.writer(open(antennafilename, "wb"), delimiter=' ')
        for ibin in range(len(antenna_size_hist)):
            antennafilewriter.writerow([ antenna_size_bins[ibin], antenna_size_hist[ibin], antenna_size_hist[ibin]/float(np.array(antenna_size_hist).sum()) ])        
        clusterfilename = self.params.output_path+"/cluster_hist_start"+str(times[0])+"_end"+str(maxtime)+"_stride"+str(stride)+".dat"
        print "writing cluster hist to", clusterfilename
        clusterfilewriter = csv.writer(open(clusterfilename, "wb"), delimiter=' ')
        for ilhcbin in range(len(cluster_1d_bins)):
            for ipsbin in range(len(cluster_1d_bins)):
                clusterfilewriter.writerow([ cluster_1d_bins[ilhcbin], cluster_1d_bins[ipsbin], cluster_2d_hist[ilhcbin][ipsbin] ])        
            clusterfilewriter.writerow([]) # write blank line so gnuplot can deal with it
            
    def PlotArrays(self):
        """public: compute and print info about PSII semi-crystalline arrays

        output file name: self.params.output_path+"/array_traj_start"+str(times[0])+"_end"+str(maxtime)+"_stride"+str(stride)+".dat"

        output file format: \n
        col1 time \n
        col2 size of largest array in any layer

        @param self The object pointer
        """
        # set up times and stride
        mintime = self.params.time
        times = self.params.getTimes()
        maxtime = np.max(times)
        stride = self.params.stride
        
        print "plotting arrays with stride", stride, "and times", times

        # set up traj file
        trajfilename = self.params.output_path+"/array_traj_start"+str(times[0])+"_end"+str(maxtime)+"_stride"+str(stride)+".dat"
        print "writing largest array trajectory to", trajfilename
        trajfilewriter = csv.writer(open(trajfilename, "wb"), delimiter=' ')

        # loop over times
        for i_time in range(len(times)):
            starttime = times[i_time]
            if starttime % stride == 0:
                print "time is", starttime
                largest_clusters = []
                # loop over layers
                for i_layer in range(len(self.trajectory[starttime][0])):
                    # build clusters, if you haven't already
                    if starttime not in self.array_sizes:
                        self.buildArrayClusters(starttime, i_layer)
                    if i_layer not in self.array_sizes[starttime]:
                        self.buildArrayClusters(starttime, i_layer)
                    #if starttime not in self._parent.data.array_sizes:
                    #    self._parent.data.array_sizes[starttime] = [ dict() for ilayer in range(len(self._parent.data.trajectory[starttime][i_type])) ]
   #                     for ilayer in range(len(self._parent.data.trajectory[starttime][i_type])):
   #                         for p in self._parent.data.trajectorystart[starttime][i_type][ilayer]:
   #                             if p:
   #                                 if p.simarrayID not in self._parent.data.array_sizes[time][ilayer]:
   #                                     self._parent.data.array_sizes[starttime][ilayer][p.simarrayID] = 1
   #                                 else:
   #                                     self._parent.data.array_sizes[starttime][ilayer][p.simarrayID] += 1
   #                         print "read", len(self._parent.data.array_sizes[starttime][ilayer].keys()), "arrays at time", starttime, "in layer", ilayer, "with largest cluster:", max(self._parent.data.array_sizes[starttime][particle._layer].values()), "and total PSII", sum(self._parent.data.array_sizes[starttime][particle._layer].values())
                    largest_clusters.append(self.array_sizes[starttime][i_layer][0][1])

                # write largest cluster to file
                trajfilewriter.writerow([starttime, max(largest_clusters)])
        print "done writing array trajectory"
        

    def DrawPS(self, time, spacing=50):
        """public: writes PostScript file with image of configuration at given time
        @param self The object pointer
        @param time Int for time
        @param spacing Number of nanometers padding to draw around each layer, default 50
        """
        # set up constants
        Lx = self.params.getWidth(time)
        Ly = self.params.getHeight(time)
        do_overlay = (2 in self.params.layers_to_draw)

        # 612x792 pixels max bounding box
        bbx = (2.0*spacing + Lx) * float(len(self.params.layers_to_draw))
        bby = 2.0*spacing + Ly
        xscale = min([1.0, bbx/612.0])
        yscale = min([1.0, bby/792.0])
        scale = min([xscale, yscale])
        bbx *= scale
        bby *= scale

        # set up output file
        sigfigs = len( str( np.max( np.array( self.trajectory.keys() ) ) ) )
        ofilename = self.params.output_path+"/frame_"
        ofilename += "layers"+''.join([str(i) for i in self.params.layers_to_draw])+"_"
        ofilename += "step"+str(time).zfill(sigfigs)+".ps" # %0*d% (time, sigfigs)

        print "writing frame to", ofilename, "with bounding box", bbx, bby
        ofile = open(ofilename, 'w')

        # print header info
        ofile.write("%!PS-Adobe-2.0 EPSF-2.0\n")
        ofile.write("%%BoundingBox: 0 0 "+str(bbx)+" "+str(bby)+"\n")
        ofile.write(str(scale)+" "+str(scale)+" scale\n")

        # abbreviated functions
        ofile.write("/tr {\n")
        ofile.write("translate\n")
        ofile.write("} def\n")

        # test if we have PSII
        if (len(self.particles[1][0]) > 0):
            sample_ps = self.particles[1][0][0]
            for layer in self.particles[1]:
                for ps in layer:
                    if ps:
                        sample_ps = ps
                        break
            if sample_ps:
                # set up PSII info
                if self.params.doEMColors:
                    prad = str(3.8)
                    plen = str(15.4)
                else:
                    prad = str(sample_ps._radius)
                    plen = str(sample_ps._length)
                xmin = str(-0.5 * sample_ps._length - sample_ps._radius)
                ymin = str(-1.0 * sample_ps._radius)
                xmax = str(0.5 * sample_ps._length + sample_ps._radius)
                ymax = str(sample_ps._radius)
                xmid = xmin
                
                # filled rod function
                ofile.write("/filled_rod{\n")
                ofile.write("newpath\n")
                ofile.write(xmid+" 0 moveto\n")
                ofile.write(xmin+" "+ymax+" "+xmax+" "+ymax+" "+prad+" arcto clear\n")
                ofile.write(xmax+" "+ymax+" "+xmax+" "+ymin+" "+prad+" arcto clear\n")
                ofile.write(xmax+" "+ymin+" "+xmin+" "+ymin+" "+prad+" arcto clear\n")
                ofile.write(xmin+" "+ymin+" "+xmin+" "+ymax+" "+prad+" arcto clear\n")
    #            ofile.write("1 0 0 setrgbcolor\n")
                ofile.write("fill\n")
                # leave out border for now
                ofile.write("} def\n")

                # draw filled rod function
                ofile.write("/dfr{\n")
                ofile.write("gsave\n")
                ofile.write("rotate filled_rod grestore\n")
                ofile.write("} def\n")

                # outline rod function
                ofile.write("/stroke_rod{\n")
                ofile.write("newpath\n")
                ofile.write(xmid+" 0 moveto\n")
                ofile.write(xmin+" "+ymax+" "+xmax+" "+ymax+" "+prad+" arcto clear\n")
                ofile.write(xmax+" "+ymax+" "+xmax+" "+ymin+" "+prad+" arcto clear\n")
                ofile.write(xmax+" "+ymin+" "+xmin+" "+ymin+" "+prad+" arcto clear\n")
                ofile.write(xmin+" "+ymin+" "+xmin+" "+ymax+" "+prad+" arcto clear\n")
    #            ofile.write("1 0 0 setrgbcolor\n")
                ofile.write("stroke\n")
                # leave out border for now
                ofile.write("} def\n")

                # draw outline rod function
                ofile.write("/dsr{\n")
                ofile.write("gsave\n")
                ofile.write("rotate stroke_rod grestore\n")
                ofile.write("} def\n")

        # test if we have LHCII
        if (len(self.particles[0][0]) > 0):
            sample_lhc = self.particles[0][0][0]
            for lhc in self.particles[0][0]:
                if lhc:
                    sample_lhc = lhc
                    break
        elif (len(self.particles[2][0])>0):
            sample_lhc = self.particles[2][0][0]
            for lhc in self.particles[2][0]:
                if lhc:
                    sample_lhc = lhc
                    break
        else:
            sample_lhc = False
        if sample_lhc:
            # set up LHCII info
            lrad = str(sample_lhc._radius)
            
            # filled circle function
            ofile.write("/circle {\n")
            ofile.write("0 0 "+lrad+" 0 360 arc\n")
            ofile.write("setrgbcolor\n")
            # leave out border for now
            ofile.write("} def\n")
            
            # abbreviated functions
            ofile.write("/cs {\n")
            ofile.write("circle stroke\n")
            ofile.write("} def\n")
            ofile.write("/cf {\n")
            ofile.write("circle fill\n")
            ofile.write("} def\n")

        # write functions for TBA-radius circles
        # filled circle function
        ofile.write("/vcf {\n")
        ofile.write("4 dict begin\n")
        ofile.write("/b exch def\n")
        ofile.write("/g exch def\n")
        ofile.write("/r exch def\n")
        ofile.write("/rad exch def\n")
        ofile.write("gsave\n")
        ofile.write("0 0 rad 0 360 arc\n")
        ofile.write("r g b setrgbcolor\n")
        ofile.write("fill\n")
        ofile.write("grestore\n")
        ofile.write("end\n")
        ofile.write("} def\n")
        # stroke circle function
        ofile.write("/vcs {\n")
        ofile.write("4 dict begin\n")
        ofile.write("/b exch def\n")
        ofile.write("/g exch def\n")
        ofile.write("/r exch def\n")
        ofile.write("/rad exch def\n")
        ofile.write("gsave\n")
        ofile.write("0 0 rad 0 360 arc\n")
        ofile.write("r g b setrgbcolor\n")
        ofile.write("stroke\n")
        ofile.write("grestore\n")
        ofile.write("end\n")
        ofile.write("} def\n")

        # fill in white background
        ofile.write("newpath\n")
        ofile.write("0 0 moveto\n")
        ofile.write("0 "+str(bby)+" lineto\n")
        ofile.write(str(bbx)+" "+str(bby)+" lineto\n")
        ofile.write(str(bbx)+" 0 lineto\n")
        ofile.write("closepath\n")
        ofile.write(str(1.0)+" setgray\n")
        ofile.write("fill\n")

        # create bounding boxes
        mid_layer_width = 0.5 * Lx
        mid_layer_height = 0.5 * Ly
        grana_rad = self.params.g_rad
        stroma_width = self.params.s_width 
        ofile.write(str(spacing)+" "+str(spacing)+" tr\n")
        for i_layer in range(len(self.params.layers_to_draw)):
            ofile.write("newpath\n")
            ofile.write(str(i_layer*(2.0*spacing + Lx))+" "+str(0)+" moveto\n")
            ofile.write(str(i_layer*(2.0*spacing + Lx))+" "+str(Ly)+" lineto\n")
            ofile.write(str(i_layer*(2.0*spacing + Lx)+Lx)+" "+str(Ly)+" lineto\n")
            ofile.write(str(i_layer*(2.0*spacing + Lx)+Lx)+" "+str(0)+" lineto\n")
            ofile.write("closepath\n")
            ofile.write("1 setlinewidth\n")
            ofile.write("0 0 0 setrgbcolor\n")
            ofile.write("stroke\n")
            ofile.write("1 setlinewidth\n")
            if stroma_width > 0 and grana_rad > 0:
                ofile.write("newpath\n")
                ofile.write(str(i_layer*(2.0*spacing + Lx))+" "+str(mid_layer_height-0.5*stroma_width)+" moveto\n")
                ofile.write(str(i_layer*(2.0*spacing + Lx))+" "+str(mid_layer_height+0.5*stroma_width)+" lineto\n")
                ofile.write(str(i_layer*(2.0*spacing + Lx)+Lx)+" "+str(mid_layer_height+0.5*stroma_width)+" lineto\n")
                ofile.write(str(i_layer*(2.0*spacing + Lx)+Lx)+" "+str(mid_layer_height-0.5*stroma_width)+" lineto\n")
                ofile.write("closepath\n")
                ofile.write("1 setlinewidth\n")
                ofile.write("0 0 0 setrgbcolor\n")
                ofile.write("stroke\n")
                ofile.write("newpath\n")
                ofile.write(str(i_layer*(2.0*spacing + Lx)+mid_layer_width-0.5*stroma_width)+" "+str(0)+" moveto\n")
                ofile.write(str(i_layer*(2.0*spacing + Lx)+mid_layer_width-0.5*stroma_width)+" "+str(Ly)+" lineto\n")
                ofile.write(str(i_layer*(2.0*spacing + Lx)+mid_layer_width+0.5*stroma_width)+" "+str(Ly)+" lineto\n")
                ofile.write(str(i_layer*(2.0*spacing + Lx)+mid_layer_width+0.5*stroma_width)+" "+str(0)+" lineto\n")
                ofile.write("closepath\n")
                ofile.write("1 setlinewidth\n")
                ofile.write("0 0 0 setrgbcolor\n")
                ofile.write("stroke\n")
                ofile.write(str(i_layer*(2.0*spacing + Lx)+mid_layer_width)+" "+str(mid_layer_height)+" "+str(grana_rad)+" 0 360 arc\n")
                ofile.write("1 1 1 setrgbcolor\n")
                ofile.write("fill\n")
                ofile.write(str(i_layer*(2.0*spacing + Lx)+mid_layer_width)+" "+str(mid_layer_height)+" "+str(grana_rad)+" 0 360 arc\n")
                ofile.write("0 0 0 setrgbcolor\n")
                ofile.write("stroke\n")

        types_to_draw = [0, 1, 2, 4, 5, 6, 7, 8]
        for i_type in types_to_draw:
            type_list = self.trajectory[time][i_type]
            for i_layer in range(len(type_list)):
                layer = type_list[i_layer]
                is_layer_to_draw = (i_layer in self.params.layers_to_draw)
                if is_layer_to_draw:
                    i_layer_for_drawing = self.params.layers_to_draw.index(i_layer)
                else:
                    i_layer_for_drawing = 0

                # set initial coords
                origin_x = i_layer_for_drawing * (Lx+2.0*spacing)
                origin_y = 0
                to_ol_x = (len(self.params.layers_to_draw) - 1 - i_layer_for_drawing) * (Lx+2.0*spacing)
                to_ol_y = 0
                ofile.write('%(x).2f %(y).2f tr\n' % {'x':origin_x, 'y':origin_y} )
                            
                for i_particle, pos in enumerate(self.trajectory[time][i_type][i_layer]):
                    if pos:
                        # get position and info
                        x = pos.x
                        y = pos.y
                        energy = pos.energy
                        color, do_draw_particle, pshape = self.getColor(i_type, energy, self.particles[i_type][i_layer][i_particle], time)
                        theta = 180.0 * pshape[0] / pi

                        if do_draw_particle and is_layer_to_draw:
                            ofile.write('%(x).2f %(y).2f tr\n' % {'x':x, 'y':y} )
                            if i_type == 1:
                                # draw PSII
                                ofile.write('%(r).2f %(g).2f %(b).2f setrgbcolor\n' %
                                            {'r':color[0], 'g':color[1], 'b':color[2]})
                                ofile.write('%(t).2f dfr\n' % {'t':theta})

#                            elif i_type in [0, 2]:
#                                # draw LHCII
#                                ofile.write('%(r).2f %(g).2f %(b).2f cf\n' %
#                                            {'r':color[0], 'g':color[1], 'b':color[2]})

                            else:
                                # draw other circle type
                                ofile.write('%(rad).2f %(r).2f %(g).2f %(b).2f vcf\n' %
                                            {'rad':pshape[1], 'r':color[0], 'g':color[1], 'b':color[2]})
                            ofile.write('%(x).2f %(y).2f tr\n' % {'x':-x, 'y':-y} )
                        
                        if do_draw_particle and do_overlay:
                            ofile.write('%(x).2f %(y).2f tr\n' % {'x':to_ol_x+x, 'y':to_ol_y+y} )
                            if i_type == 1:
                                # draw PSII
                                ofile.write('%(r).2f %(g).2f %(b).2f setrgbcolor\n' %
                                            {'r':color[0], 'g':color[1], 'b':color[2]})
                                ofile.write('%(t).2f dsr\n' % {'t':theta})
#                            elif i_type in [0, 2]:
#                                # draw LHCII
#                                ofile.write('%(r).2f %(g).2f %(b).2f cs\n' %
#                                            {'r':color[0], 'g':color[1], 'b':color[2]})
                            else:
                                # draw other circle type
                                ofile.write('%(rad).2f %(r).2f %(g).2f %(b).2f vcs\n' %
                                            {'rad':pshape[1], 'r':color[0], 'g':color[1], 'b':color[2]})
                            ofile.write('%(x).2f %(y).2f tr\n' % {'x':-to_ol_x-x, 'y':-to_ol_y-y} )
                ofile.write('%(x).2f %(y).2f tr\n' % {'x':-origin_x, 'y':-origin_y} )
                
        # finish up
        ofile.write("showpage\n")
        ofile.close()
        print "done writing frame"


    def getColor(self, i_type, energy, particle, time):
        """public: get RGB color for particle based on parameters in self.params
        @param self The object pointer
        @param i_type Int for particle type
        @param energy particle energy
        @param particle Particle object
        @param time Int for time
        @retval tuple RGB values between 0 and 1
        @retval bool True if particle will be drawn, False if not
        @retval tuple (theta, pradius, plength), which are modified for PSII in doEMColors
        """
        
        # set up return values
        color = (0, 0, 0)
        do_draw_particle = True
        stride = self.params.stride
        pos = self.trajectory[time][i_type][particle._layer][particle._id]
        pshape = (pos.theta, particle._radius, particle._length)

        if i_type == 0 or i_type == 2: # LHCII, free or s.c.
            if self.params.doEMColors:
                do_draw_particle = False
            elif self.params.doTrimerEnergyColors:
                frac = 1 - abs(energy / self.params.lhc_stacking_epsilon)
                if frac < 0 or frac > 1:
                    print "epsilon is wrong, got energy", energy
                else:
                    color = self.params.getInterpolatedColor(i_type, frac)

            elif self.params.doDispColors:
                if i_type == 0:
                    disp = 0
                    maxdisp = 2.0 * 2.0*particle._radius
                    if time+stride in self.params.times:
                        stridepos = self.trajectory[time+stride][i_type][particle._layer][particle._id]
                        currpos = self.trajectory[time][i_type][particle._layer][particle._id]
                        if stridepos:
                            dispx = self.PBC_diff(stridepos.x, currpos.x, self.params.getWidth(time))
                            dispy = self.PBC_diff(stridepos.y, currpos.y, self.params.getHeight(time))
                            disp = sqrt(dispx*dispx + dispy*dispy)
                            if disp > maxdisp:
                                disp = maxdisp
                            #if time not in particle.displacements:
                            #    particle.displacements[time] = disp
                    color = colors.colorConverter.to_rgb(cm.jet(disp/maxdisp))
                else:
                    do_draw_particle = False

            elif self.params.doClusterColors:
                if i_type == 0:
                    if not self.trajectory[time][i_type][particle._layer][particle._id].already_clustered:
                        self.buildETClusters(time, particle._layer)
                    icluster = self.trajectory[time][i_type][particle._layer][particle._id].clusterID
                    colorval = float(self.cluster_sizes[time][particle._layer][icluster][1])
                    if colorval > 0:
                        colorval = (self.cluster_sizes[time][particle._layer][icluster][0] + 2.0*self.cluster_sizes[time][particle._layer][icluster][1])/float(self.cluster_sizes[time][particle._layer][icluster][1])
                        color = colors.colorConverter.to_rgb(cm.jet(colorval/10.0))                       
                    elif len(self.trajectory[time][1][particle._layer]) > 0:
                        color = colors.colorConverter.to_rgb(cm.jet(0))                       
                    else:
                        color = colors.colorConverter.to_rgb(cm.jet(float(icluster) / float(len(self.cluster_sizes[time][particle._layer]))))    
                else:
                    do_draw_particle = False

            elif self.params.doMTrimerColors:
                do_draw_particle = True
                if i_type == 0:
                    if particle._tagged:
                        color = self.params.MTrimerColor
                    else:
                        color = self.params.freeTrimerColor
                else:
                    color = self.params.STrimerColor

            elif self.params.doXtalColors: # do crystal colors
                if i_type == 0:
                    color = self.params.lhcDefaultColor
                else:
                    do_draw_particle = False
                    
            elif self.params.doHexaticColors:
                hexatic = self.calculateHexatic(i_type, particle._id, particle._layer, time)
                color = colors.colorConverter.to_rgb(cm.jet(hexatic))
       
            elif self.params.doTagColors:
                if particle._tagged:
                    color = self.params.lhcTagColor
                else:
                    color = self.params.lhcDefaultColor

            elif self.params.doPhosColors and particle.isPhos():
                color = self.params.lhcTagColor

            else:
                color = self.params.lhcDefaultColor
        elif i_type == 1: # PSII
            if self.params.doEMColors:
              #  theta -= (2*i_layer - 1) * 0.389
              #  plength = 15.3
              #  pradius = 3.8
              #  pshape = (pos.theta - (2*particle._layer - 1) * 0.389, 15.3, 3.8)
                color = (1.0, 1.0, 1.0)
            elif self.params.doPSStackColors:
                frac = 1 - abs(energy / self.params.ps_stacking_epsilon)
                if frac < 0 or frac > 1:
                    print "epsilon is wrong, got energy", energy
                    color = (0,0,0)
                else:
                    color = self.params.getInterpolatedColor(i_type, frac)
            elif particle._tagged and self.params.doTagColors:
                color = self.params.psTagColor
            elif self.params.doXtalColors:
                if not self.trajectory[time][i_type][particle._layer][particle._id].already_arrayed:
                    self.buildArrayClusters(time, particle._layer)
                icluster = self.trajectory[time][i_type][particle._layer][particle._id].arrayID
             #   colorval = log(float(self.array_sizes[time][particle._layer][icluster][1])) / log(float(self.array_sizes[time][particle._layer][0][1]))
             #   color = colors.colorConverter.to_rgb(cm.jet(colorval))
                colorval = self.array_sizes[time][particle._layer][icluster][1]
                if colorval > 10:
                    color = self.params.psXtalColor
                else:
                    color = self.params.psFluidColor
            elif self.params.doDelaunayColors:
                colorval = self.calulateDelaunayShapes(i_type, particle._id, particle._layer, time)
                color = colors.colorConverter.to_rgb(cm.jet(colorval))    
            elif self.params.doDispColors:
                    disp = 0
                    maxdisp = 2.0 * 2.0*particle._radius
                    if time+stride in self.params.times:
                        stridepos = self.trajectory[time+stride][i_type][particle._layer][particle._id]
                        currpos = self.trajectory[time][i_type][particle._layer][particle._id]
                        if stridepos:
                            dispx = self.PBC_diff(stridepos.x, currpos.x, self.params.getWidth(time))
                            dispy = self.PBC_diff(stridepos.y, currpos.y, self.params.getHeight(time))
                            disp = sqrt(dispx*dispx + dispy*dispy)
                            if disp > maxdisp:
                                disp = maxdisp
                            #if time not in particle.displacements:
                            #    particle.displacements[time] = disp
                    color = colors.colorConverter.to_rgb(cm.jet(disp/maxdisp))
            elif self.params.doClusterColors:
                if not self.trajectory[time][i_type][particle._layer][particle._id].already_clustered:
                     self.buildETClusters(time, particle._layer)
                icluster = self.trajectory[time][i_type][particle._layer][particle._id].clusterID
                colorval = float(self.cluster_sizes[time][particle._layer][icluster][1])
                if colorval > 0: colorval = (self.cluster_sizes[time][particle._layer][icluster][0] + 2.0*self.cluster_sizes[time][particle._layer][icluster][1])/float(self.cluster_sizes[time][particle._layer][icluster][1])
                color = colors.colorConverter.to_rgb(cm.jet(colorval/10.0))
            else:
                color = self.params.psDefaultColor
        elif i_type == 3: # LHCII monomer
            do_draw_particle = False
        elif i_type == 4: # megacomplex site
            if self.params.doMsiteEnergyColors:
                frac = 1 - abs(energy / self.params.megacomplex_epsilon)
                if frac < 0 or frac > 1:
                    print "epsilon is wrong, got energy", energy
                    color = (0,0,0)
                else:
                    color = self.params.getInterpolatedColor(i_type, frac)
            else:
                do_draw_particle = False
        elif i_type in [5, 6, 7, 8]: # generic hard disc
            if self.params.doDispColors:
                disp = 0
                maxdisp = 2.0 * 2.0*particle._radius
                if time+stride in self.params.times:
                    stridepos = self.trajectory[time+stride][i_type][particle._layer][particle._id]
                    currpos = self.trajectory[time][i_type][particle._layer][particle._id]
                    if stridepos:
                        dispx = self.PBC_diff(stridepos.x, currpos.x, self.params.getWidth(time))
                        dispy = self.PBC_diff(stridepos.y, currpos.y, self.params.getHeight(time))
                        disp = sqrt(dispx*dispx + dispy*dispy)
                        if disp > maxdisp:
                            disp = maxdisp
                        #if time not in particle.displacements:
                        #    particle.displacements[time] = disp
                color = colors.colorConverter.to_rgb(cm.jet(disp/maxdisp))
            elif self.params.doTagColors and particle._tagged:
                color = self.params.discTagColor
            elif self.params.doDelaunayColors:
                hexatic = self.calculateDelaulayColors(i_type, particle._id, particle._layer, time)
                color = colors.colorConverter.to_rgb(cm.jet(hexatic))
              #  if hexatic > 0.5:
              #      color = self.params.discDefaultColor
              #  else:
              #      color = self.params.discTagColor
            else:
                color = self.params.discDefaultColor

        return color, do_draw_particle, pshape

class Datapoint:
    """Class for storing single particle, single time information"""
    def __init__(self, x=0.0, y=0.0, theta=0.0, energy=0.0, region=0, cluster=-1):
        """public: constructor
        @param self The object pointer
        """
        # initialized from file
        self.x = x
        self.y = y
        self.theta = theta
        self.energy = energy
        self.region = region
        self.simarrayID = cluster

        # will compute later if wanted
        self.nbhrs = []
        self.nematic = False
        self.hexatic = False
        self.end1_x = False
        self.end1_y = False
        self.end2_x = False
        self.end2_y = False
        self.already_clustered = False
        self.clusterID = -1
        self.already_arrayed = False
        self.arrayID = -1
        self.in_linear_array = False

    def findEnds(self, halflength):
        if not self.end1_x:
            self.end1_x = halflength * cos(self.theta) + self.x
        if not self.end1_y:
            self.end1_y = halflength * sin(self.theta) + self.y
        if not self.end2_x:
            self.end2_x = halflength * cos(self.theta+pi) + self.x
        if not self.end2_y:
            self.end2_y = halflength * sin(self.theta+pi) + self.y

    def update_images(self, images, lastpos, Lx, Ly):
        # update x image
        if self.x < 0.1 * Lx and lastpos.x > 0.9 * Lx:
            images[0] = images[0] + 1
        elif self.x > 0.9 * Lx and lastpos.x < 0.1 * Lx:
            images[0] = images[0] - 1

        # update y image
        if self.y < 0.1 * Ly and lastpos.y > 0.9 * Ly:
            images[1] = images[1] + 1
        elif self.y > 0.9 * Ly and lastpos.y < 0.1 * Ly:
            images[1] = images[1] - 1

        return images

class Particle:
    def __init__(self, type, id, layer, radius, length, phos):
        self._type = type
        self._id = id
        self._layer = layer
        self._radius = radius
        self._length = length
        self._phos = phos
        self._tagged = False
        self.displacements = dict()
        
    def getID(self):
        return (self._type, self._layer, self._id)
        
    def isPhos(self):
        return self._phos

class Params:
    def __init__(self):
        # initialize with default parameters
        self.width = 300
        self.height = 300
        self.g_rad = 0
        self.s_width = 0
        self.time = 0
        self.times = [0]
        self.stride = 100000
        self.doPBC = True
        self.output_path = "./"
        self.system_sizes = dict()
        self.system_sizes[0] = [self.width, self.height, self.g_rad, self.s_width] 
        self.layers_to_draw = [0, 1, 2]

        self.lhc_stacking_epsilon = 3.0
        self.ps_stacking_epsilon = 0.0
        self.megacomplex_epsilon = 2.0

        # target crystal type
        self.crystal = "c2s2mboekema"

        # colors from colorbrewer2.org, qualitative paired
        self.lhcDefaultColor = (178.0/255.0, 223.0/255.0, 138.0/255.0)
        self.lhcTagColor = (51.0/255.0, 160.0/255.0, 44.0/255.0)
        self.psDefaultColor = (116.0/255.0, 206.0/255.0, 227.0/255.0)
        self.psTagColor = (31.0/255.0, 120.0/255.0, 180.0/255.0)
        self.discDefaultColor = (253.0/255.0, 191.0/255.0, 111.0/255.0)
        self.discTagColor = (255.0/255.0, 127.0/255.0, 0.0/255.0)

        # colors from colorbrewer2.org, sequential green 6-color        
        self.freeTrimerColor = (161.0/255.0, 217.0/255.0, 155.0/255.0)
        self.MTrimerColor = (49.0/255.0, 163.0/255.0, 84.0/255.0)
        self.STrimerColor = (0.0/255.0, 109.0/255.0, 44.0/255.0)

        # colors from colorbrewer2.org, qualitative set1
        self.psXtalColor = (228.0/255.0, 26.0/255.0, 28.0/255.0)
        self.psFluidColor = (55.0/255.0, 126.0/255.0, 184.0/255.0)
        
        # my LHC energy scheme
        self.lhcEnergyMin = (0.0/255.0, 100.0/255.0, 0.0/255.0)
        self.lhcEnergyMax = (220.0/255.0, 150.0/255.0, 0.0/255.0)
        self.lhcPhos = (220.0/255.0, 75.0/255.0, 0.0/255.0)
        self.psEnergyMin = (0.0/255.0, 100.0/255.0, 0.0/255.0)
        self.psEnergyMax = (220.0/255.0, 150.0/255.0, 0.0/255.0)
        self.megaEnergyMin = (255.0/255.0, 127.0/255.0, 0.0/255.0)
        self.megaEnergyMax = (255.0/255.0, 255.0/255.0, 255.0/255.0)
 
        # color booleans
        self.doTrimerEnergyColors = True
        self.doMonomerEnergyColors = False
        self.doMsiteEnergyColors = False
        self.doMTrimerColors = False
        self.doTagColors = True
        self.doPSStackColors = False
        self.doXtalColors = False
        self.doDelaunayColors = False
        self.doEMColors = False
        self.doDispColors = False
        self.doClusterColors = False
        self.doHexaticColors = False
        self.doPhosColors = True

    def addTime(self, newtime):
        if newtime in self.times:
            print newtime, "already in times:", self.times
        else:
            self.times.append(newtime)
            self.times.sort()

    def setTime(self, newtime):
      #  if newtime not in self.times:
      #      self.addTime(newtime)
        self.time = newtime
        self.setWidth(self.system_sizes[newtime][0])
        self.setHeight(self.system_sizes[newtime][1])
        self.g_rad = self.system_sizes[newtime][2]
        self.s_width = self.system_sizes[newtime][3]

    def getTimes(self):
        nptimes = np.array(self.times)
        return nptimes[ np.where( nptimes >= self.time ) ]

    def addSystemSize(self, newtime, newwidth, newheight, grana_rad=0, stroma_width=0):
        self.system_sizes[newtime] = [newwidth, newheight, grana_rad, stroma_width]

    def setWidth(self, newwidth):
        self.width = newwidth
        if self.width != self.system_sizes[self.time][0]:
            print "width is set to", self.width, "but value grabbed from filename was", self.system_sizes[self.time][0]
            self.system_sizes[self.time][0] = self.width

    def getWidth(self, t):
        return self.system_sizes[t][0]

    def setHeight(self, newheight):
        self.height = newheight
        if self.height != self.system_sizes[self.time][1]:
            print "height is set to", self.height, "but value grabbed from filename is",self.system_sizes[self.time][1]
            self.system_sizes[self.time][1] = self.height

    def getHeight(self, t):
        return self.system_sizes[t][1]

    def setGRad(self, newgrad):
        self.g_rad = newgrad
        if self.g_rad != self.system_sizes[self.time][2]:
            print "grana radius is set to", self.g_rad, "but value grabbed from filename is",self.system_sizes[self.time][2]
            self.system_sizes[self.time][2] = self.g_rad

    def getGRad(self, t):
        return self.system_sizes[t][2]

    def setSWidth(self, newswidth):
        self.s_width = newswidth
        if self.s_width != self.system_sizes[self.time][3]:
            print "stroma width is set to", self.s_width, "but value grabbed from filename is",self.system_sizes[self.time][3]
            self.system_sizes[self.time][3] = self.s_width

    def getSWidth(self, t):
        return self.system_sizes[t][3]

    def setStride(self, newstride):
        self.stride = newstride

    def setStrideAuto(self):
        try:
            self.stride = np.min(np.diff(self.times))
        except ValueError:
            self.stride = 0

    def setPBC(self, pbcbool):
        self.doPBC = pbcbool

    def setOutputPath(self, pathstring):
        self.output_path = pathstring

    def getInterpolatedColor(self, type, fraction):
        if type == 0 or type == 2:
            return (self.lhcEnergyMin[0] + fraction*(self.lhcEnergyMax[0]-self.lhcEnergyMin[0]),
                    self.lhcEnergyMin[1] + fraction*(self.lhcEnergyMax[1]-self.lhcEnergyMin[1]),
                    self.lhcEnergyMin[2] + fraction*(self.lhcEnergyMax[2]-self.lhcEnergyMin[2])
                    )
        elif type == 1:
            return (self.psEnergyMin[0] + fraction*(self.psEnergyMax[0]-self.psEnergyMin[0]),
                    self.psEnergyMin[1] + fraction*(self.psEnergyMax[1]-self.psEnergyMin[1]),
                    self.psEnergyMin[2] + fraction*(self.psEnergyMax[2]-self.psEnergyMin[2])
                    )
        elif type == 4:
            return (self.megaEnergyMin[0] + fraction*(self.megaEnergyMax[0]-self.megaEnergyMin[0]),
                    self.megaEnergyMin[1] + fraction*(self.megaEnergyMax[1]-self.megaEnergyMin[1]),
                    self.megaEnergyMin[2] + fraction*(self.megaEnergyMax[2]-self.megaEnergyMin[2])
                    )
        else:
            print "can only interpolate for types 0,1,2,4"
            return (0,0,0)
