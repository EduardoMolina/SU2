#!/usr/bin/env python

## \file 
#  \brief
#  \author Eduardo Molina
#  \version 7.0.8 "Blackbird"
#
# The current SU2 release has been coordinated by the
# SU2 International Developers Society <www.su2devsociety.org>
# with selected contributions from the open-source community.
#
# The main research teams contributing to the current release are:
#  - Prof. Juan J. Alonso's group at Stanford University.
#  - Prof. Piero Colonna's group at Delft University of Technology.
#  - Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
#  - Prof. Alberto Guardone's group at Polytechnic University of Milan.
#  - Prof. Rafael Palacios' group at Imperial College London.
#  - Prof. Vincent Terrapon's group at the University of Liege.
#  - Prof. Edwin van der Weide's group at the University of Twente.
#  - Lab. of New Concepts in Aeronautics at Tech. Institute of Aeronautics.
#
# Copyright 2012-2020, Francisco D. Palacios, Thomas D. Economon,
#                      Tim Albring, and the SU2 contributors.
#
# SU2 is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
#
# SU2 is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with SU2. If not, see <http://www.gnu.org/licenses/>.

# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------

import sys
from optparse import OptionParser # use a parser for configuration
import pysu2                  # imports the SU2 wrapped module
from pybemt.solver import Solver
#import SUAVE
from math import *
import numpy as np
from scipy.interpolate import interp1d
import pdb

# -------------------------------------------------------------------
#  Read Actuator Disk File 
# -------------------------------------------------------------------

def ReadActuatorDiskFile(filename):

  infile = open(filename, 'r')
  lines  = infile.readlines()
  infile.close()

  up_marker   = []
  down_marker = []
  ad_center   = []
  ad_axis     = []
  ad_radius   = []
  ad_J        = []

  for line in lines:
    if "MARKER_ACTDISK=" in line:
      up_marker.append(line.strip().split()[1])
      down_marker.append(line.strip().split()[2])
    if "CENTER=" in line:
      aux = line.strip().split()[1:]
      ad_center.append([float(i) for i in aux])      
    if "AXIS=" in line:
      aux = line.strip().split()[1:]
      ad_axis.append([float(i) for i in aux])
    if "RADIUS=" in line:
      aux = line.strip().split()[1]
      ad_radius.append(float(aux))
    if "ADV_RATIO=" in line:
      aux = line.strip().split()[1]
      ad_J.append(float(aux))


  ActDiskDict   = {
    "Up_Marker"   : up_marker,
    "Down_Marker" : down_marker,
    "nProp"       : len(up_marker),
    "Up_ID"       : [-1] * len(up_marker),
    "Up_nVertex"  : [-1] * len(up_marker),
    "Down_ID"     : [-1] * len(up_marker),
    "Down_nVertex": [-1] * len(up_marker),
    "Center"      : ad_center,
    "Axis"        : ad_axis,
    "Radius"      : ad_radius,
    "J"           : ad_J
  } #specified by the user

  return ActDiskDict

def Interpolate(ActDiskDict, Coord, ):

   return
# -------------------------------------------------------------------
#  Main 
# -------------------------------------------------------------------

def main():

  # Command line options
  parser=OptionParser()
  parser.add_option("-f", "--file", dest="filename", help="Read config from FILE", metavar="FILE")
  parser.add_option("--nDim", dest="nDim", default=2, help="Define the number of DIMENSIONS",
                    metavar="DIMENSIONS")
  parser.add_option("--nZone", dest="nZone", default=1, help="Define the number of ZONES",
                    metavar="ZONES")
  parser.add_option("--parallel", action="store_true",
                    help="Specify if we need to initialize MPI", dest="with_MPI", default=False)

  (options, args) = parser.parse_args()
  options.nDim  = int( options.nDim )
  options.nZone = int( options.nZone )

  # Import mpi4py for parallel run
  has_mpi = False
  if options.with_MPI == True:
    from mpi4py import MPI
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    has_mpi = True
  else:
    comm = 0
    rank = 0

  if (rank == 0):
    print("\nPython Wrapper BEMT\n")

  # Initialize the corresponding driver of SU2, this includes solver preprocessing
  try:
    SU2Driver = pysu2.CSinglezoneDriver(options.filename, options.nZone, comm);
  except TypeError as exception:
    print('A TypeError occured in pysu2.CDriver : ',exception)
    if options.with_MPI == True:
      print('ERROR : You are trying to initialize MPI with a serial build of the wrapper. Please, remove the --parallel option that is incompatible with a serial build.')
    else:
      print('ERROR : You are trying to launch a computation without initializing MPI but the wrapper has been built in parallel. Please add the --parallel option in order to initialize MPI for the wrapper.')
    return


  # Get all the markers defined on this rank and their associated indices.
  allMarkerIDs = SU2Driver.GetAllBoundaryMarkers()
  print(allMarkerIDs) 
  # Read actuator disk file and create the dictionary
  ActDiskDict = ReadActuatorDiskFile("ActuatorDisk.dat")
  nProp = ActDiskDict["nProp"]

  # Instantiate the external BEMT solver
  bemt = Solver("prowim_optimized.ini")

  # Check if the specified marker exists and the number of vertices in this rank.
  for item in allMarkerIDs.items():
    for i in range(nProp):
      if ActDiskDict["Up_Marker"][i] == item[0]:
        ActDiskDict["Up_ID"][i] = item[1]
        ActDiskDict["Up_nVertex"][i] = SU2Driver.GetNumberVertices(item[1])
      if ActDiskDict["Down_Marker"][i] == item[0]:
        ActDiskDict["Down_ID"][i] = item[1]
        ActDiskDict["Down_nVertex"][i] = SU2Driver.GetNumberVertices(item[1])
  #print(ActDiskDict)

  # Retrieve some control parameters from the driver
  TimeIter  = SU2Driver.GetTime_Iter()
  nTimeIter = SU2Driver.GetnInner_Iter()
  #nTimeIter = 1000
  nBEMT     = 20

  # Time loop is defined in Python so that we have acces to SU2 functionalities at each time step
  if rank == 0:
    print("\n------------------------------ Begin Solver -----------------------------\n")
  sys.stdout.flush()

  if options.with_MPI == True:
    comm.Barrier()

  if rank == 0:
    print("\n------------------------------ Reading ActDisk File -----------------------------\n")
  SU2Driver.ReadActuatordDiskFile()

  if options.with_MPI == True:
    comm.Barrier()
    
  # Time iteration preprocessing
  SU2Driver.Preprocess(TimeIter)
  
  for iBEMT in range(1,nBEMT):

    if rank == 0:
      print("\nExternal BEMT iteration", iBEMT)

    # Get actuator disk integrated data from SU2_CFD
    for i in range(nProp):

      VelActDisk = SU2Driver.GetActDiskInflowVelocity(ActDiskDict["Up_Marker"][i])
      RhoActDisk = SU2Driver.GetActDiskInflowDensity(ActDiskDict["Up_Marker"][i])
      ViscActDisk = SU2Driver.GetActDiskInflowViscosity(ActDiskDict["Up_Marker"][i])

      #print('Rank ', rank, VelActDisk, RhoActDisk, ViscActDisk)
      # Set the free-stream velocity, density and viscosity on the external BEMT solver.
      bemt.fluid.rho = RhoActDisk
      bemt.fluid.mu  = ViscActDisk
      bemt.v_inf     = VelActDisk
        
      # Run the BEMT solver and get dCt/dr_ and dCP/dr_ (note that r_ = r/R)
      T_bemt, Q_bemt, P_bemt, df_bemt    = bemt.run()
      J_bemt, CT_bemt, CQ_bemt, CP_bemt, eta_bemt = bemt.rotor_coeffs(T_bemt, Q_bemt, P_bemt)
      if (rank == 0):
        print('--- pyBEMT Results ---')
        print('Advance Ratio :\t',J_bemt)
        print('Thrust (N):\t',T_bemt)
        print('Torque (Nm):\t',Q_bemt)
        print('Power (W):\t',P_bemt)
        print('CT :\t',CT_bemt)
        print('CQ :\t',CQ_bemt)
        print('CP :\t',CP_bemt)
        print('eta :\t',eta_bemt)

      R_, dCT_, dCP_ = bemt.force_per_unit_area()
        
      Fa = dCT_*(2.*RhoActDisk*pow(VelActDisk,2))/(pow(J_bemt,2)*np.pi*R_)
      Ft = dCP_*(2*RhoActDisk*pow(VelActDisk,2))/((J_bemt*np.pi*R_) * (J_bemt*np.pi*R_))

      # R_ = np.linspace(0.2031, 0.9795, 37)
      # Fa = np.ones_like(R_)
      # Ft = np.ones_like(R_)

      # Build the linear interpolator
      f_Ct = interp1d(R_, Fa, fill_value=(Fa[0],Fa[-1]), bounds_error= False)
      f_Cp = interp1d(R_, Ft, fill_value=(Ft[0],Ft[-1]), bounds_error= False)

      if ActDiskDict["Up_nVertex"][i] > 0:

        # Interpolate to actuator disk variable load (Inlet AD)
        for j in range(ActDiskDict["Up_nVertex"][i]):
         
          ad_center = np.array(ActDiskDict["Center"][i])
          ad_radius = ActDiskDict["Radius"][i]

          # Get the coordinates of the current node
          xp       = SU2Driver.GetVertexCoordX(ActDiskDict["Up_ID"][i], j)
          yp       = SU2Driver.GetVertexCoordY(ActDiskDict["Up_ID"][i], j)
          zp       = SU2Driver.GetVertexCoordZ(ActDiskDict["Up_ID"][i], j)
          P        = np.array([xp, yp, zp])

          # Computation of the radius coordinates for the current node
          P_radius  = P - ad_center

          #Computation of the non-dimensional radius for the current node.
          P_radius_ = np.linalg.norm(P_radius) / ad_radius

          #Fx, Fy and Fz are the x, y and z components of the tangential and radial forces per unit area resultant. 
          Fa_i = 0.0; Fx_i = 0.0; Fy_i = 0.0; Fz_i = 0.0;
          if (P_radius_ > 0.):
            Fa_i = f_Ct(P_radius_)
            Ft_i = f_Cp(P_radius_)

            Fx_i = Ft_i * P_radius[0] / (P_radius_ * ad_radius)
            Fy_i = Ft_i * P_radius[2] / (P_radius_ * ad_radius)
            Fz_i = -Ft_i * P_radius[1] / (P_radius_ * ad_radius)
          
          SU2Driver.SetVertexActuatorDiskForce(ActDiskDict["Up_ID"][i], j, 0, float(Fa_i))
          SU2Driver.SetVertexActuatorDiskForce(ActDiskDict["Up_ID"][i], j, 1, float(Fx_i))
          SU2Driver.SetVertexActuatorDiskForce(ActDiskDict["Up_ID"][i], j, 2, float(Fy_i))
          SU2Driver.SetVertexActuatorDiskForce(ActDiskDict["Up_ID"][i], j, 3, float(Fz_i))

      if ActDiskDict["Down_nVertex"][i] > 0:
        # Interpolate to actuator disk variable load (Outlet AD)
        for j in range(ActDiskDict["Down_nVertex"][i]):
         
          ad_center = np.array(ActDiskDict["Center"][i])
          ad_radius = ActDiskDict["Radius"][i]

          # Get the coordinates of the current node
          xp       = SU2Driver.GetVertexCoordX(ActDiskDict["Down_ID"][i], j)
          yp       = SU2Driver.GetVertexCoordY(ActDiskDict["Down_ID"][i], j)
          zp       = SU2Driver.GetVertexCoordZ(ActDiskDict["Down_ID"][i], j)
          P        = np.array([xp, yp, zp])

          # Computation of the radius coordinates for the current node
          P_radius  = P - ad_center

          #Computation of the non-dimensional radius for the current node.
          P_radius_ = np.linalg.norm(P_radius) / ad_radius

          #Fx, Fy and Fz are the x, y and z components of the tangential and radial forces per unit area resultant. 
          Fa_i = 0.0; Fx_i = 0.0; Fy_i = 0.0; Fz_i = 0.0;
          if (P_radius_ > 0.):
            Fa_i = f_Ct(P_radius_)
            Ft_i = f_Cp(P_radius_)

            Fx_i = Ft_i * P_radius[0] / (P_radius_ * ad_radius)
            Fy_i = Ft_i * P_radius[2] / (P_radius_ * ad_radius)
            Fz_i = -Ft_i * P_radius[1] / (P_radius_ * ad_radius)
          
          SU2Driver.SetVertexActuatorDiskForce(ActDiskDict["Down_ID"][i], j, 0, float(Fa_i))
          SU2Driver.SetVertexActuatorDiskForce(ActDiskDict["Down_ID"][i], j, 1, float(Fx_i))
          SU2Driver.SetVertexActuatorDiskForce(ActDiskDict["Down_ID"][i], j, 2, float(Fy_i))
          SU2Driver.SetVertexActuatorDiskForce(ActDiskDict["Down_ID"][i], j, 3, float(Fz_i))
      
    # Run the simulation
    SU2Driver.Run()
    
    TimeIter = iBEMT * nTimeIter

    # Update the solver for the next time iteration
    #SU2Driver.Update()

    # Monitor the solver and output solution to file if required
    stopCalc = SU2Driver.Monitor(TimeIter)
    SU2Driver.Output(TimeIter)
    #if (stopCalc == True):
    #  break

  # Postprocess the solver and exit cleanly
  SU2Driver.Postprocessing()

  if SU2Driver != None:
    del SU2Driver

# -------------------------------------------------------------------
#  Run Main Program
# -------------------------------------------------------------------

# this is only accessed if running from command prompt
if __name__ == '__main__':
    main()
