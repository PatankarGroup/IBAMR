c
c     Routines to compute discrete stresses on patches.
c
c     Created on 21 June 2016 by Nishant Nangia
c
c     Copyright (c) 2002-2014, Boyce Griffith
c     All rights reserved.
c
c     Redistribution and use in source and binary forms, with or without
c     modification, are permitted provided that the following conditions
c     are met:
c
c        * Redistributions of source code must retain the above
c          copyright notice, this list of conditions and the following
c          disclaimer.
c
c        * Redistributions in binary form must reproduce the above
c          copyright notice, this list of conditions and the following
c          disclaimer in the documentation and/or other materials
c          provided with the distribution.
c
c        * Neither the name of The University of North Carolina nor the
c          names of its contributors may be used to endorse or promote
c          products derived from this software without specific prior
c          written permission.
c
c     THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND
c     CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,
c     INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
c     MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
c     DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS
c     BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
c     EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
c     TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
c     DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
c     ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR
c     TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF
c     THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
c     SUCH DAMAGE.
c
define(NDIM,3)dnl
define(REAL,`double precision')dnl
define(INTEGER,`integer')dnl
include(SAMRAI_FORTDIR/pdat_m4arrdim3d.i)dnl


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Computes w = stress u.
c
c     Uses centered differences to compute the edge centered stress of a
c     side centered vector field u=(u0,u1,u2).
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine stoestress3d(
     &     W0,W1,W2,W_gcw,
     &     U0,U1,U2,U_gcw,
     &     ilower0,iupper0,
     &     ilower1,iupper1,
     &     ilower2,iupper2,
     &     dx)
c
      implicit none
c
c     Input.
c
      INTEGER ilower0,iupper0
      INTEGER ilower1,iupper1
      INTEGER ilower2,iupper2
      INTEGER W_gcw,U_gcw

      REAL U0(SIDE3d0(ilower,iupper,U_gcw))
      REAL U1(SIDE3d1(ilower,iupper,U_gcw))
      REAL U2(SIDE3d2(ilower,iupper,U_gcw))

      REAL dx(0:NDIM-1)
c
c     Input/Output.
c
      REAL W0(EDGE3d0(ilower,iupper,W_gcw))
      REAL W1(EDGE3d1(ilower,iupper,W_gcw))
      REAL W2(EDGE3d2(ilower,iupper,W_gcw))
c
c     Local variables.
c
      INTEGER i0,i1,i2
c
c     Compute the edge centered stress of u=(u0,u1,u2).
c

      do i2 = ilower2,iupper2+1
         do i1 = ilower1,iupper1+1
            do i0 = ilower0,iupper0
               W0(i0,i1,i2) = (U2(i0,i1,i2) - U2(i0,i1-1,i2)  )/dx(1)  
     &                      + (U1(i0,i1,i2) - U1(i0,i1  ,i2-1))/dx(2)
     
               W0(i0,i1,i2) = 0.5*W0(i0,i1,i2)
            enddo
         enddo
      enddo

      do i2 = ilower2,iupper2+1
         do i1 = ilower1,iupper1
            do i0 = ilower0,iupper0+1
               W1(i0,i1,i2) = (U0(i0,i1,i2) - U0(i0,i1,i2-1))/dx(2)  
     &                      + (U2(i0,i1,i2) - U2(i0-1,i1,i2))/dx(0)
               W1(i0,i1,i2) = 0.5*W1(i0,i1,i2)
            enddo
         enddo
      enddo

      do i2 = ilower2,iupper2
         do i1 = ilower1,iupper1+1
            do i0 = ilower0,iupper0+1
               W2(i0,i1,i2) = (U1(i0,i1,i2) - U1(i0-1,i1,i2))/dx(0) 
     &                      + (U0(i0,i1,i2) - U0(i0,i1-1,i2))/dx(1)
               W2(i0,i1,i2) = 0.5*W2(i0,i1,i2)
            enddo
         enddo
      enddo
c
      return
      end
c