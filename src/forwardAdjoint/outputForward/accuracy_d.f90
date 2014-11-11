   !        Generated by TAPENADE     (INRIA, Tropics team)
   !  Tapenade 3.10 (r5363) -  9 Sep 2014 09:53
   !
   !
   !      ******************************************************************
   !      *                                                                *
   !      * File:          inputParam.f90                                  *
   !      * Author:        Edwin van der Weide, Steve Repsher,             *
   !      *                C.A.(Sandy) Mader                               *
   !      * Starting date: 12-11-2002                                      *
   !      * Last modified: 09-17-2009                                      *
   !      *                                                                *
   !      ******************************************************************
   !
   MODULE ACCURACY_D
   !
   !      ******************************************************************
   !      *                                                                *
   !      * Definition of some parameters which make the code more         *
   !      * readable. The actual values of this parameters are arbitrary;  *
   !      * in the code always the symbolic names are (should be) used.    *
   !      *                                                                *
   !      ******************************************************************
   !
   USE PRECISION
   IMPLICIT NONE
   SAVE 
   !
   INTEGER(kind=inttype), PARAMETER :: firstorder=1, secondorder=2, &
   & thirdorder=3, fourthorder=4, fifthorder=5
   END MODULE ACCURACY_D
