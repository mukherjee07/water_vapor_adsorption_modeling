"""
Author: Krishnendu Mukherjee, Chemical Engineering, University of Notre Dame
Advisor: Dr. Yamil Colon
Goal: The objective of this code is to create a CIF and a block file (i.e. crystallographic information file) which has the following features:
NOTE: Read this section before using this code. 
Input --> C-C bond length (They would be cylinder made from Carbon element), Number of outer ring (1,2,3 in the radial direction), Number of carbon layers (in the z-direction), Number of Carbon atoms in one single ring
Block File -->
Features: It would have 2 number of blockades
1st type -- square above and below the cylinder; thus this would be implemented on a rectangular section; positions
of x - [0,0.2,0.4,0.6,0.8,1.0], y - [0,0.2,0.4,0.6,0.8,1.0], z - { [1+z2]/2, z1/2 };z1 and z2 would need to be extracted
for implementing this
2nd type -- cylinderical; from z1 to z2. x - [xo + (r + atom_size/2)*sin(theta) + del] where theta --> (0,2 pi in 20 steps)
y would follow the same protocol but xo would be replaced by yo, and sine by cosine
z --> zi + c*(z2 - z1)/N where N == z1 - z2/bond_length, and c --> [0,N]
and thus we expect to generate both the .cif and block (no extension) using this code.
The code is designed in such a way that it can be used with varying number of carbon rings, cell length and others
However, if the system is not exclusively carbon or perhaps not cylinderical; many things needs to be changed.
Infact a new script would be required.

"""

#importing libraries
import numpy as np
import pandas as pd
import math
import os

#defining key input variables
bond_length = 1.54 #In Angstroem
layers = 10
rings = 2
atomsperring = XX
filename = "IC_YY_AA"
blockfilename = "IC_YY_AA.block"
C = 0.7
delta = 0.5
inter_ring_distance = 1.54

#defining the cell length parameters with respect to layers variable
cell_length = 40
z_length = 2*bond_length*layers

#defining flag variable and the origin position
layer_counter = 0
counter = 0
xo = cell_length/2 # @center
yo = cell_length/2 # @center
zo = bond_length/2 # beginning from bottom of the surfacei; just enough so the carbon atom touches the bottom wall
z_position = zo
layer_counter = 0

#creating the xyz file using shell cmd
os.system("rm -f "+str(filename)+".cif")
os.system("touch "+str(filename)+".cif")

#adding the header of a general cif file type and upadating the cell length
os.system("cat header.cif >> "+str(filename)+".cif")
os.system("sed -i 's/PARAM/"+str(cell_length)+"/' "+str(filename)+".cif ")
os.system("sed -i 's/ZLENGTH/"+str(z_length)+"/' "+str(filename)+".cif ")

#Generating coordinates in x,y, and z fashion and then converting them to fractional coordinates
while z_position < z_length:
   for ring_counter in range(rings):
      r = (atomsperring*(bond_length/2)/math.pi) + ring_counter*inter_ring_distance
      for position_counter in range(atomsperring):
         counter += 1
         theta = 360*(position_counter/atomsperring)
         # Generating coordinates in cartesian coordinates
         temp_x = xo + r*math.sin(math.radians(theta))
         temp_y = yo + r*math.cos(math.radians(theta))
         temp_z = zo + (layer_counter*bond_length)
         # Now converting to fractional coordinates
         temp_x = temp_x/cell_length
         temp_y = temp_y/cell_length
         temp_z = temp_z/z_length
         if ring_counter == 0:
            if (layer_counter % 2) == 0:
               charge = 0.3915
            else:
               charge = -0.3915
         else:
            if (layer_counter % 2) == 0:
               charge = -0.3915
            else:
               charge = 0.3915
         #print("C",temp_x,temp_y,temp_z)
         print(temp_z)
         if temp_x < 1 and temp_y < 1 and temp_z < 1:
            os.system("echo "+"C"+str(counter)+" "+"C "+str(temp_x)+" "+str(temp_y)+" "+str(temp_z)+" "+str(charge)+" >> "+str(filename)+".cif")
   layer_counter += 1
   z_position += bond_length


#Creating a block file to prevent any active simulation on vacant spaces
os.system("touch "+str(blockfilename))

#Populating the block file with spherical 'blocks' to fill the void space in the cif
#Dealing with square blocks; (easier perhaps :P)

#checking if the blockfile is present or not; if there old one would be removed and a new instance would be created
os.system("rm -f "+str(blockfilename))
os.system("echo TotalBlocks >> "+str(blockfilename))

#Creating important variables
z1 = 0
z2 = 1
count_s = 0 #counter for number of square blocks
N_grid = 12 #counter for number of grids on one axis

'''
#starting the loop for square blocks
for x in range(N_grid):
   for y in range(N_grid):
      for z in range(2):
         if z == 1:
            os.system("echo "+str(x/N_grid)+" "+str(y/N_grid)+" "+str((1+z2)/2)+" "+str(((1-z2)/2)*cell_length)+" >> "+str(blockfilename))
            count_s += 1
         if z == 0:
            os.system("echo "+str(x/N_grid)+" "+str(y/N_grid)+" "+str((z1)/2)+" "+str((z1/2)*cell_length)+" >> "+str(blockfilename))
            count_s += 1
'''
##Dealing with cylinderical blocks, more tricky ;()

#Defining key parameters
N = (z_length)/C
N = math.floor(N)
print(N)
count_c = 1 
z_position = 0
#Xe_dia = 

# NOTE: If you want to test with Ur atom then simply uncomment the two lines with Ur on it
#starting the loop for cylindrical blocks
#Note: this loop moves radially outwards until it hits the wall (the 'if' statement would take care of that)

for z_n in range(N):
   for ring_c in range(60):
      for theta in range(72):
         if ring_c ==0:
            x_position = (xo+((r+C+delta)*math.sin(math.radians(theta*5))))/cell_length
            y_position = (yo+((r+C+delta)*math.cos(math.radians(theta*5))))/cell_length
            z_position =  (z_n*C)/z_length
            if (x_position < 1 and y_position < 1 and x_position > 0 and y_position > 0) and (((x_position-xo)**2 + (y_position-yo)**2) > ((r+C)**2)): 
               os.system("echo "+str(x_position)+" "+str(y_position)+" "+str(z_position)+" "+str(C)+" >> "+str(blockfilename))
               #os.system("echo "+"Xe"+str(count_c)+" "+"Xe "+str(x_position)+" "+str(y_position)+" "+str(z_position)+" >> "+str(filename)+".cif")
               count_c += 1
               print(count_c)
         elif ring_c == 1 or ring_c == 2 or ring_c ==3:
            x_position = (xo+(((ring_c+1)*C+r)*math.sin(math.radians(theta*5))))/cell_length
            y_position = (yo+(((ring_c+1)*C+r)*math.cos(math.radians(theta*5))))/cell_length
            z_position = (z_n*C)/z_length
            if (x_position < 1 and y_position < 1 and x_position > 0 and y_position > 0) and (((x_position-xo)**2 + (y_position-yo)**2) > ((r+((3+ring_c)*C)))**2): 
               os.system("echo "+str(x_position)+" "+str(y_position)+" "+str(z_position)+" "+str(2*C)+" >> "+str(blockfilename))
               #Block with Uranium testing               
               #os.system("echo "+"Xe"+str(count_c)+" "+"Xe "+str(x_position)+" "+str(y_position)+" "+str(z_position)+" >> "+str(filename)+".cif")
               count_c += 1
         else:
            x_position = (xo+(((ring_c+1)*C+r)*math.sin(math.radians(theta*5))))/cell_length
            y_position = (yo+(((ring_c+1)*C+r)*math.cos(math.radians(theta*5))))/cell_length
            z_position = (z_n*C)/z_length
            if (x_position < 1 and y_position < 1 and x_position > 0 and y_position > 0) and (((x_position-xo)**2 + (y_position-yo)**2) > ((r+((3+ring_c)*C)))**2): 
               os.system("echo "+str(x_position)+" "+str(y_position)+" "+str(z_position)+" "+str(3*C)+" >> "+str(blockfilename))
               #Block with Uranium testing               
               #os.system("echo "+"Xe"+str(count_c)+" "+"Xe "+str(x_position)+" "+str(y_position)+" "+str(z_position)+" >> "+str(filename)+".cif")
               count_c += 1
	
#adding blocks in the 8 corners

os.system("echo 1.0 0.0 0.875 "+str(5*C)+" >> "+str(blockfilename))
os.system("echo 1.0 0.0 0.875 "+str(5*C)+" >> "+str(blockfilename))
os.system("echo 0.0 0.0 0.875 "+str(5*C)+" >> "+str(blockfilename))
os.system("echo 0.0 1.0 0.875 "+str(5*C)+" >> "+str(blockfilename))
os.system("echo 1.0 1.0 0.625 "+str(5*C)+" >> "+str(blockfilename))
os.system("echo 0.0 0.0 0.625 "+str(5*C)+" >> "+str(blockfilename))
os.system("echo 0.0 1.0 0.625 "+str(5*C)+" >> "+str(blockfilename))
os.system("echo 1.0 1.0 0.625 "+str(5*C)+" >> "+str(blockfilename))
os.system("echo 0.0 0.0 0.750 "+str(5*C)+" >> "+str(blockfilename))
os.system("echo 0.0 1.0 0.750 "+str(5*C)+" >> "+str(blockfilename))
os.system("echo 1.0 0.0 0.750 "+str(5*C)+" >> "+str(blockfilename))
os.system("echo 1.0 1.0 0.750 "+str(5*C)+" >> "+str(blockfilename))
os.system("echo 1.0 0.0 0.125 "+str(5*C)+" >> "+str(blockfilename))
os.system("echo 0.0 0.0 0.125 "+str(5*C)+" >> "+str(blockfilename))
os.system("echo 0.0 1.0 0.125 "+str(5*C)+" >> "+str(blockfilename))
os.system("echo 1.0 1.0 0.125 "+str(5*C)+" >> "+str(blockfilename))
os.system("echo 0.0 0.0 0.3750 "+str(5*C)+" >> "+str(blockfilename))
os.system("echo 0.0 1.0 0.3750 "+str(5*C)+" >> "+str(blockfilename))
os.system("echo 1.0 0.0 0.3750 "+str(5*C)+" >> "+str(blockfilename))
os.system("echo 1.0 1.0 0.3750 "+str(5*C)+" >> "+str(blockfilename))
os.system("echo 0.0 0.0 0.25 "+str(5*C)+" >> "+str(blockfilename))
os.system("echo 0.0 1.0 0.25 "+str(5*C)+" >> "+str(blockfilename))
os.system("echo 1.0 0.0 0.25 "+str(5*C)+" >> "+str(blockfilename))
os.system("echo 1.0 1.0 0.25 "+str(5*C)+" >> "+str(blockfilename))
os.system("echo 0.0 0.0 0.50 "+str(5*C)+" >> "+str(blockfilename))
os.system("echo 0.0 1.0 0.50 "+str(5*C)+" >> "+str(blockfilename))
os.system("echo 1.0 0.0 0.50 "+str(5*C)+" >> "+str(blockfilename))
os.system("echo 1.0 1.0 0.50 "+str(5*C)+" >> "+str(blockfilename))
os.system("echo 0.0 0.0 0.0 "+str(5*C)+" >> "+str(blockfilename))
os.system("echo 0.0 1.0 0.0 "+str(5*C)+" >> "+str(blockfilename))
os.system("echo 0.0 0.0 1.0 "+str(5*C)+" >> "+str(blockfilename))
os.system("echo 0.0 1.0 1.0 "+str(5*C)+" >> "+str(blockfilename))
os.system("echo 1.0 0.0 0.0 "+str(5*C)+" >> "+str(blockfilename))
os.system("echo 1.0 0.0 1.0 "+str(5*C)+" >> "+str(blockfilename))
os.system("echo 1.0 1.0 0.0 "+str(5*C)+" >> "+str(blockfilename))
os.system("echo 1.0 1.0 1.0 "+str(5*C)+" >> "+str(blockfilename))
count_c += 36

#using 'sed' command in shell to replace number of blocks with actual number of total blocks
total_blocks = count_s + count_c
os.system("sed -i 's/TotalBlocks/"+str(total_blocks)+"/' "+str(blockfilename))

#DONE

