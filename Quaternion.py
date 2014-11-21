# -*- coding: ascii -*-

"""
Quaternion library with numpy support
(adapted from Alex Holkner's pyeuclid library)
~~~~~~~~~~~~~

#
# This library is free software; you can redistribute it and/or modify it
# under the terms of the GNU Lesser General Public License as published by the
# Free Software Foundation; either version 2.1 of the License, or (at your
# option) any later version.
# 
# This library is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
# FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License
# for more details.
# 
# You should have received a copy of the GNU Lesser General Public License
# along with this library; if not, write to the Free Software Foundation,
# Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301 USA

"""
import numpy as np

class Quaternion:
    __slots__ = ['w', 'x', 'y', 'z']

	def __init__(self, w=np.array([1]), x=np.array([0]), y=np.array([0]), z=np.array([0])):
        self.w = w
        self.x = x
        self.y = y
        self.z = z

    def __copy__(self):
        Q = Quaternion()
        Q.w = self.w
        Q.x = self.x
        Q.y = self.y
        Q.z = self.z
        return Q
		
    def euler(self):
        R=np.empty((3,3,max(self.w.shape)))
        R[0,0,:] = 2*np.square(self.w)-1+2*np.square(self.x);
        R[1,0,:] = 2*( np.multiply(self.x,self.y)- np.multiply(self.w,self.z));
        R[2,0,:] = 2*( np.multiply(self.x,self.z)+ np.multiply(self.w,self.y));
        R[2,1,:] = 2*( np.multiply(self.y,self.z)- np.multiply(self.w,self.x));
        R[2,2,:] = 2*np.square(self.w)-1+2*np.square(self.z);
        phi = np.arctan2(R[2,1,:], R[2,2,:] );
        theta = -np.arctan( np.divide(R[2,0,:],np.sqrt(1-np.square(R[2,0,:]))) );
        psi = np.arctan2(R[1,0,:], R[0,0,:] );
        eul=np.column_stack((phi[:].T,theta[:].T,psi[:].T))
        return eul
		
    copy = __copy__

    def __abs__(self):
        return np.sqrt(np.square(self.w)+np.square(self.x)+np.square(self.y)+np.square(self.z))
		
	magnitude=__abs__
	
    def __mul__(self, other):
        if isinstance(other, Quaternion):
            Ax = self.x
            Ay = self.y
            Az = self.z
            Aw = self.w
            Bx = other.x
            By = other.y
            Bz = other.z
            Bw = other.w
            Q = Quaternion()
            Q.x =  np.multiply(Ax,Bw) + np.multiply(Ay,Bz) - np.multiply(Az,By) + np.multiply(Aw,Bx)    
            Q.y = np.multiply(-Ax,Bz) + np.multiply(Ay,Bw) + np.multiply(Az,Bx) + np.multiply(Aw,By)
            Q.z =  np.multiply(Ax,By) - np.multiply(Ay,Bx) + np.multiply(Az,Bw) + np.multiply(Aw,Bz)
            Q.w = np.multiply(-Ax,Bx) - np.multiply(Ay,By) - np.multiply(Az,Bz) + np.multiply(Aw,Bw)
            return Q	
			
    def conjugate(self):
        self.x = -self.x
        self.y = -self.y
        self.z = -self.z

    def conjugated(self):
        Q=Quaternion(w=self.w,x=self.x,y=self.y,z=self.z)
		Q.conjugate()
		return Q

    def normalize(self):
        d = self.magnitude()
        self.w = np.divide(self.w,d)
        self.x = np.divide(self.x,d)
        self.y = np.divide(self.y,d)
        self.z = np.divide(self.z,d)

    def normalized(self):
        Q = Quaternion(w=self.w,x=self.x,y=self.y,z=self.z)
		Q.normalize()
		return Q
